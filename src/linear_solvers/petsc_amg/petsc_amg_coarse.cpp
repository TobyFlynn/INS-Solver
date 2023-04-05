#include "linear_solvers/petsc_amg_coarse.h"

#include "linear_solvers/petsc_utils.h"
#include "timing.h"

#include <iostream>
#include <type_traits>

extern Timing *timer;

PETScAMGCoarseSolver::PETScAMGCoarseSolver(DGMesh *m) : PETScAMGSolver(m) {

}

bool PETScAMGCoarseSolver::solve(op_dat rhs, op_dat ans) {
  timer->startTimer("PETScAMGCoarseSolver - solve");
  timer->startTimer("PETScAMGCoarseSolver - get PETSc matrix");
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
  timer->endTimer("PETScAMGCoarseSolver - get PETSc matrix");

  if(bc)
    matrix->apply_bc(rhs, bc);

  Vec b, x;
  PETScUtils::create_vec_p_adapt(&b, matrix->getUnknowns());
  PETScUtils::create_vec_p_adapt(&x, matrix->getUnknowns());

  PETScUtils::load_vec_coarse(&b, rhs);
  PETScUtils::load_vec_coarse(&x, ans);

  timer->startTimer("PETScAMGCoarseSolver - KSPSolve");
  KSPSolve(ksp, b, x);
  timer->endTimer("PETScAMGCoarseSolver - KSPSolve");

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
  PETScUtils::store_vec_coarse(&solution, ans);

  PETScUtils::destroy_vec(&b);
  PETScUtils::destroy_vec(&x);

  timer->endTimer("PETScAMGCoarseSolver - solve");

  return converged;
}
