#include "linear_solvers/petsc_jacobi.h"

#ifdef INS_MPI
#include "mpi_helper_func.h"
#include <iostream>
#include "op_mpi_core.h"
#endif

PetscErrorCode matAMultJacobi(Mat A, Vec x, Vec y) {
  PETScJacobiSolver *poisson;
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

// TODO update for p-adaptivity
void PETScJacobiSolver::create_shell_mat() {
  if(pMatInit)
    MatDestroy(&pMat);

  // const int mat_size = matrix->getUnknowns();
  const int mat_size = mesh->cells->size * DG_NP;
  MatCreateShell(PETSC_COMM_WORLD, mat_size, mat_size, PETSC_DETERMINE, PETSC_DETERMINE, this, &pMat);
  MatShellSetOperation(pMat, MATOP_MULT, (void(*)(void))matAMultJacobi);
  MatShellSetVecType(pMat, VECSTANDARD);

  pMatInit = true;
}

PetscErrorCode preconJacobi(PC pc, Vec x, Vec y) {
  PETScJacobiSolver *poisson;
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

void PETScJacobiSolver::set_shell_pc(PC pc) {
  PCShellSetApply(pc, preconJacobi);
  PCShellSetContext(pc, this);
}
