#include "linear_solvers/petsc_pmultigrid.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

PetscErrorCode matAMultPM(Mat A, Vec x, Vec y) {
  PETScPMultigrid *poisson;
  MatShellGetContext(A, &poisson);
  const DG_FP *x_ptr;
  DG_FP *y_ptr;
  VecGetArrayRead(x, &x_ptr);
  VecGetArray(y, &y_ptr);

  poisson->calc_rhs(x_ptr, y_ptr);

  VecRestoreArrayRead(x, &x_ptr);
  VecRestoreArray(y, &y_ptr);
  return 0;
}

void PETScPMultigrid::create_shell_mat() {
  if(pMatInit)
    MatDestroy(&pMat);

  MatCreateShell(PETSC_COMM_WORLD, matrix->getUnknowns(), matrix->getUnknowns(), PETSC_DETERMINE, PETSC_DETERMINE, this, &pMat);
  MatShellSetOperation(pMat, MATOP_MULT, (void(*)(void))matAMultPM);
  MatShellSetVecType(pMat, VECSTANDARD);

  pMatInit = true;
}

PetscErrorCode preconPM(PC pc, Vec x, Vec y) {
  PETScPMultigrid *poisson;
  PCShellGetContext(pc, (void **)&poisson);
  const DG_FP *x_ptr;
  DG_FP *y_ptr;
  VecGetArrayRead(x, &x_ptr);
  VecGetArray(y, &y_ptr);

  poisson->precond(x_ptr, y_ptr);

  VecRestoreArrayRead(x, &x_ptr);
  VecRestoreArray(y, &y_ptr);
  return 0;
}

void PETScPMultigrid::set_shell_pc(PC pc) {
  PCShellSetApply(pc, preconPM);
  PCShellSetContext(pc, this);
}
