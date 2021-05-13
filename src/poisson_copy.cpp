#include "poisson.h"

// Copy u PETSc vec array to OP2 dat (TODO avoid this copy)
void Poisson_MF2::copy_u(const double *u_d) {
  op_arg u_copy_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 1, u_copy_args);
  memcpy(u->data, u_d, u->set->size * 15 * sizeof(double));
  op_mpi_set_dirtybit(1, u_copy_args);
}

// Copy rhs OP2 dat to PETSc vec array (TODO avoid this copy)
void Poisson_MF2::copy_rhs(double *rhs_d) {
  op_arg rhs_copy_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges(data->cells, 1, rhs_copy_args);
  memcpy(rhs_d, rhs->data, rhs->set->size * 15 * sizeof(double));
  op_mpi_set_dirtybit(1, rhs_copy_args);
}

// Create a PETSc vector for CPUs
void Poisson::create_vec(Vec *v, int size) {
  // VecCreateSeq(PETSC_COMM_SELF, size * data->cells->size, v);
  VecCreate(PETSC_COMM_WORLD, v);
  VecSetType(*v, VECSTANDARD);
  VecSetSizes(*v, size * data->cells->size, PETSC_DECIDE);
}

// Destroy a PETSc vector
void Poisson::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for CPUs
void Poisson::load_vec(Vec *v, op_dat v_dat, int size) {
  double *v_ptr;
  VecGetArray(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, size, "double", OP_READ)
  };
  op_mpi_halo_exchanges(data->cells, 1, vec_petsc_args);
  memcpy(v_ptr, (double *)v_dat->data, size * v_dat->set->size * sizeof(double));
  op_mpi_set_dirtybit(1, vec_petsc_args);
  VecRestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for CPUs
void Poisson::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecGetArrayRead(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 1, vec_petsc_args);
  memcpy((double *)v_dat->data, v_ptr, 15 * v_dat->set->size * sizeof(double));
  op_mpi_set_dirtybit(1, vec_petsc_args);
  VecRestoreArrayRead(*v, &v_ptr);
}

// Create a PETSc matrix for CPUs
void Poisson::create_mat(Mat *m, int row, int col, int prealloc) {
  // MatCreate(PETSC_COMM_SELF, m);
  MatCreate(PETSC_COMM_WORLD, m);
  // MatSetSizes(*m, PETSC_DECIDE, PETSC_DECIDE, row, col);
  MatSetSizes(*m, row, col, PETSC_DECIDE, PETSC_DECIDE);
  // MatSetType(*m, MATSEQAIJ);
  MatSetType(*m, MATAIJ);
  #ifdef INS_MPI
  MatMPIAIJSetPreallocation(*m, prealloc, NULL, prealloc, NULL);
  #else
  MatSeqAIJSetPreallocation(*m, prealloc, NULL);
  #endif
}

PetscErrorCode matAMult2(Mat A, Vec x, Vec y) {
  timer->startLinearSolveMFMatMult();
  Poisson_MF2 *poisson;
  MatShellGetContext(A, &poisson);
  const double *x_ptr;
  double *y_ptr;
  VecGetArrayRead(x, &x_ptr);
  VecGetArray(y, &y_ptr);

  poisson->calc_rhs(x_ptr, y_ptr);

  VecRestoreArrayRead(x, &x_ptr);
  VecRestoreArray(y, &y_ptr);
  timer->endLinearSolveMFMatMult();
  return 0;
}

void Poisson_MF2::create_shell_mat(Mat *m) {
  // MatCreateShell(PETSC_COMM_SELF, 15 * data->cells->size, 15 * data->cells->size, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatCreateShell(PETSC_COMM_WORLD, 15 * data->cells->size, 15 * data->cells->size, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatShellSetOperation(*m, MATOP_MULT, (void(*)(void))matAMult2);
  MatShellSetVecType(*m, VECSTANDARD);
}
