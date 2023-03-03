#include "linear_solvers/petsc_block_jacobi.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  PETScBlockJacobiSolver *poisson;
  MatShellGetContext(A, &poisson);
  const DG_FP *x_ptr;
  DG_FP *y_ptr;
  VecCUDAGetArrayRead(x, &x_ptr);
  VecCUDAGetArray(y, &y_ptr);

  poisson->calc_rhs(x_ptr, y_ptr);

  VecCUDARestoreArrayRead(x, &x_ptr);
  VecCUDARestoreArray(y, &y_ptr);
  return 0;
}
// TODO update for p-adaptivity
void PETScBlockJacobiSolver::create_shell_mat() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreateShell(PETSC_COMM_WORLD, matrix->getUnknowns(), matrix->getUnknowns(), PETSC_DETERMINE, PETSC_DETERMINE, this, &pMat);
  MatShellSetOperation(pMat, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(pMat, VECCUDA);

  pMatInit = true;
}

PetscErrorCode precon(PC pc, Vec x, Vec y) {
  PETScBlockJacobiSolver *poisson;
  PCShellGetContext(pc, (void **)&poisson);
  const DG_FP *x_ptr;
  DG_FP *y_ptr;
  VecCUDAGetArrayRead(x, &x_ptr);
  VecCUDAGetArray(y, &y_ptr);

  poisson->precond(x_ptr, y_ptr);

  VecCUDARestoreArrayRead(x, &x_ptr);
  VecCUDARestoreArray(y, &y_ptr);
  return 0;
}

void PETScBlockJacobiSolver::set_shell_pc(PC pc) {
  PCShellSetApply(pc, precon);
  PCShellSetContext(pc, this);
}
