#include "linear_solvers/petsc_inv_mass.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#endif

PetscErrorCode matAMultInvMass(Mat A, Vec x, Vec y) {
  PETScInvMassSolver *poisson;
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
void PETScInvMassSolver::create_shell_mat() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreateShell(PETSC_COMM_WORLD, matrix->getUnknowns(), matrix->getUnknowns(), PETSC_DETERMINE, PETSC_DETERMINE, this, &pMat);
  MatShellSetOperation(pMat, MATOP_MULT, (void(*)(void))matAMultInvMass);
  MatShellSetVecType(pMat, VECCUDA);

  pMatInit = true;
}

PetscErrorCode preconInvMass(PC pc, Vec x, Vec y) {
  PETScInvMassSolver *poisson;
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

void PETScInvMassSolver::set_shell_pc(PC pc) {
  PCShellSetApply(pc, preconInvMass);
  PCShellSetContext(pc, this);
}
