#include "poisson.h"

// Copy u PETSc vec array to OP2 dat (TODO avoid this copy)
void Poisson_MF::copy_u(const double *u_d) {
  op_arg u_copy_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, u_copy_args);
  cudaMemcpy(u->data_d, u_d, data->numCells * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, u_copy_args);
}

// Copy rhs OP2 dat to PETSc vec array (TODO avoid this copy)
void Poisson_MF::copy_rhs(double *rhs_d) {
  op_arg rhs_copy_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, rhs_copy_args);
  cudaMemcpy(rhs_d, rhs->data_d, data->numCells * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, rhs_copy_args);
}

// Copy u PETSc vec array to OP2 dat (TODO avoid this copy)
void Poisson_MF2::copy_u(const double *u_d) {
  op_arg u_copy_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, u_copy_args);
  cudaMemcpy(u->data_d, u_d, data->numCells * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, u_copy_args);
}

// Copy rhs OP2 dat to PETSc vec array (TODO avoid this copy)
void Poisson_MF2::copy_rhs(double *rhs_d) {
  op_arg rhs_copy_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, rhs_copy_args);
  cudaMemcpy(rhs_d, rhs->data_d, data->numCells * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, rhs_copy_args);
}

// Create a PETSc vector for GPUs
void Poisson::create_vec(Vec *v, int size) {
  VecCreateSeqCUDA(PETSC_COMM_SELF, size * data->numCells, v);
}

// Destroy a PETSc vector
void Poisson::destroy_vec(Vec *v) {
  VecDestroy(v);
}

// Load a PETSc vector with values from an OP2 dat for GPUs
void Poisson::load_vec(Vec *v, op_dat v_dat, int size) {
  double *v_ptr;
  VecCUDAGetArray(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, size, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, vec_petsc_args);
  cudaMemcpy(v_ptr, (double *)v_dat->data_d, size * data->numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, vec_petsc_args);
  VecCUDARestoreArray(*v, &v_ptr);
}

// Load an OP2 dat with the values from a PETSc vector for GPUs
void Poisson::store_vec(Vec *v, op_dat v_dat) {
  const double *v_ptr;
  VecCUDAGetArrayRead(*v, &v_ptr);
  op_arg vec_petsc_args[] = {
    op_arg_dat(v_dat, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, vec_petsc_args);
  cudaMemcpy((double *)v_dat->data_d, v_ptr, 15 * data->numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, vec_petsc_args);
  VecCUDARestoreArrayRead(*v, &v_ptr);
}

// Create a PETSc matrix for GPUs
void Poisson::create_mat(Mat *m, int row, int col, int prealloc) {
  MatCreate(PETSC_COMM_SELF, m);
  MatSetSizes(*m, PETSC_DECIDE, PETSC_DECIDE, row, col);
  MatSetType(*m, MATSEQAIJCUSPARSE);
  MatSeqAIJSetPreallocation(*m, prealloc, NULL);
}

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  timer->startLinearSolveMFMatMult();
  Poisson_MF *poisson;
  MatShellGetContext(A, &poisson);
  const double *x_ptr;
  double *y_ptr;
  VecCUDAGetArrayRead(x, &x_ptr);
  VecCUDAGetArray(y, &y_ptr);

  poisson->calc_rhs(x_ptr, y_ptr);

  VecCUDARestoreArrayRead(x, &x_ptr);
  VecCUDARestoreArray(y, &y_ptr);
  timer->endLinearSolveMFMatMult();
  return 0;
}

PetscErrorCode matAMult2(Mat A, Vec x, Vec y) {
  timer->startLinearSolveMFMatMult();
  Poisson_MF2 *poisson;
  MatShellGetContext(A, &poisson);
  const double *x_ptr;
  double *y_ptr;
  VecCUDAGetArrayRead(x, &x_ptr);
  VecCUDAGetArray(y, &y_ptr);

  poisson->calc_rhs(x_ptr, y_ptr);

  VecCUDARestoreArrayRead(x, &x_ptr);
  VecCUDARestoreArray(y, &y_ptr);
  timer->endLinearSolveMFMatMult();
  return 0;
}

void Poisson_MF::create_shell_mat(Mat *m) {
  MatCreateShell(PETSC_COMM_SELF, 15 * data->numCells, 15 * data->numCells, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatShellSetOperation(*m, MATOP_MULT, (void(*)(void))matAMult);
  MatShellSetVecType(*m, VECCUDA);
}

void Poisson_MF2::create_shell_mat(Mat *m) {
  MatCreateShell(PETSC_COMM_SELF, 15 * data->numCells, 15 * data->numCells, PETSC_DETERMINE, PETSC_DETERMINE, this, m);
  MatShellSetOperation(*m, MATOP_MULT, (void(*)(void))matAMult2);
  MatShellSetVecType(*m, VECCUDA);
}
