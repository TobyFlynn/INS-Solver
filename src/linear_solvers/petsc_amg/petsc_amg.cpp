#include "linear_solvers/petsc_amg.h"

#include "linear_solvers/petsc_utils.h"

#include <iostream>

PETScAMGSolver::PETScAMGSolver(DGMesh2D *m) {
  mesh = m;
  nullspace = false;
  pMatInit = false;

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  KSPSetTolerances(ksp, 1e-10, 1e-50, 1e5, 2.5e2);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCGAMG);
  PCGAMGSetNSmooths(pc, 4);
  PCGAMGSetSquareGraph(pc, 1);
  PCGAMGSetNlevels(pc, 20);
  PCMGSetLevels(pc, 20, NULL);
  PCMGSetCycleType(pc, PC_MG_CYCLE_W);
  PCGAMGSetRepartition(pc, PETSC_TRUE);
  PCGAMGSetReuseInterpolation(pc, PETSC_TRUE);
}

PETScAMGSolver::~PETScAMGSolver() {
  KSPDestroy(&ksp);
}

bool PETScAMGSolver::solve(op_dat rhs, op_dat ans) {
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

  if(bc)
    matrix->apply_bc(rhs, bc);

  Vec b, x;
  PETScUtils::create_vec_p_adapt(&b, matrix->unknowns);
  PETScUtils::create_vec_p_adapt(&x, matrix->unknowns);

  PETScUtils::load_vec_p_adapt(&b, rhs, mesh);
  PETScUtils::load_vec_p_adapt(&x, ans, mesh);

  KSPSolve(ksp, b, x);

  int numIt;
  KSPGetIterationNumber(ksp, &numIt);
  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  // Check that the solver converged
  bool converged = true;
  if(reason < 0) {
    double residual;
    KSPGetResidualNorm(ksp, &residual);
    converged = false;
    std::cout << "Number of iterations for linear solver: " << numIt << std::endl;
    std::cout << "Converged reason: " << reason << " Residual: " << residual << std::endl;
  }

  Vec solution;
  KSPGetSolution(ksp, &solution);
  PETScUtils::store_vec_p_adapt(&solution, ans, mesh);

  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);

  return converged;
}
