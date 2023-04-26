#include "linear_solvers/petsc_amg.h"

#include "linear_solvers/petsc_utils.h"
#include "timing.h"
#include "config.h"

#include <iostream>
#include <type_traits>

extern Timing *timer;
extern Config *config;

PETScAMGSolver::PETScAMGSolver(DGMesh *m) {
  bc = nullptr;
  mesh = m;
  nullspace = false;
  pMatInit = false;

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  // KSPSetType(ksp, KSPGMRES);
  KSPSetType(ksp, KSPCG);
  if(std::is_same<DG_FP,double>::value)
    KSPSetTolerances(ksp, 1e-8, 1e-50, 1e5, 2.5e2);
  else
    KSPSetTolerances(ksp, 1e-5, 1e-50, 1e5, 2.5e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCGAMG);
  // PCGAMGSetNSmooths(pc, 4);
  // PCGAMGSetSquareGraph(pc, 1);
  int tmp = 20;
  config->getInt("petsc-amg", "levels", tmp);
  PCGAMGSetNlevels(pc, tmp);
  PCMGSetLevels(pc, tmp, NULL);
  tmp = 1;
  config->getInt("petsc-amg", "w-cycle", tmp);
  PCMGSetCycleType(pc, tmp == 1 ? PC_MG_CYCLE_W : PC_MG_CYCLE_V);
  tmp = 1;
  config->getInt("petsc-amg", "repartition", tmp);
  PCGAMGSetRepartition(pc, tmp == 1 ? PETSC_TRUE : PETSC_FALSE);
  tmp = 1;
  config->getInt("petsc-amg", "reuse-interpolation", tmp);
  PCGAMGSetReuseInterpolation(pc, tmp == 1 ? PETSC_TRUE : PETSC_FALSE);
  tmp = 0;
  config->getInt("petsc-amg", "coarse-eq-lim", tmp);
  if(tmp > 0) {
    PCGAMGSetCoarseEqLim(pc, tmp);
  }
  tmp = 0;
  config->getInt("petsc-amg", "proc-eq-lim", tmp);
  if(tmp > 0) {
    PCGAMGSetProcEqLim(pc, tmp);
  }
  tmp = 0;
  config->getInt("petsc-amg", "use-cpu-for-coarse-solve", tmp);
  PCGAMGSetCpuPinCoarseGrids(pc, tmp == 1 ? PETSC_TRUE : PETSC_FALSE);
  double tmp_d = 0.0;
  config->getDouble("petsc-amg", "threshold-scale", tmp_d);
  if(tmp_d > 0) {
    PCGAMGSetThresholdScale(pc, tmp_d);
  }
  tmp = 0;
  config->getInt("petsc-amg", "use-aggs", tmp);
  PCGAMGASMSetUseAggs(pc, tmp == 1 ? PETSC_TRUE : PETSC_FALSE);
}

PETScAMGSolver::~PETScAMGSolver() {
  KSPDestroy(&ksp);
}

bool PETScAMGSolver::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("PETScAMGSolver - solve");
  timer->startTimer("PETScAMGSolver - get PETSc matrix");
  if(matrix->getPETScMat(&pMat)) {
    if(nullspace) {
      MatNullSpace ns;
      MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &ns);
      MatSetNullSpace(*pMat, ns);
      MatSetTransposeNullSpace(*pMat, ns);
      MatNullSpaceDestroy(&ns);
    }
    KSPSetOperators(ksp, *pMat, *pMat);
  }
  timer->endTimer("PETScAMGSolver - get PETSc matrix");

  if(bc)
    matrix->apply_bc(rhs, bc);

  Vec b, x;
  PETScUtils::create_vec_p_adapt(&b, matrix->getUnknowns());
  PETScUtils::create_vec_p_adapt(&x, matrix->getUnknowns());

  PETScUtils::load_vec_p_adapt(&b, rhs, mesh);
  PETScUtils::load_vec_p_adapt(&x, ans, mesh);

  timer->startTimer("PETScAMGSolver - KSPSolve");
  KSPSolve(ksp, b, x);
  timer->endTimer("PETScAMGSolver - KSPSolve");

  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  // Check that the solver converged
  bool converged = true;
  if(reason < 0) {
    DG_FP residual;
    KSPGetResidualNorm(ksp, &residual);
    converged = false;
    std::cout << "Number of iterations for AMG linear solver: " << numIt << std::endl;
    std::cout << "Converged reason: " << reason << " Residual: " << residual << std::endl;
  }

  Vec solution;
  KSPGetSolution(ksp, &solution);
  PETScUtils::store_vec_p_adapt(&solution, ans, mesh);

  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);

  timer->endTimer("PETScAMGSolver - solve");

  return converged;
}

void PETScAMGSolver::set_tol(const DG_FP tol) {
  KSPSetTolerances(ksp, tol, 1e-50, 1e5, 2.5e2);
}
