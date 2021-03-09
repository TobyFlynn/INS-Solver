#include "poisson.h"

void Poisson::copy_u(const double *u) {
  op_arg u_copy_args[] = {
    op_arg_dat(pU, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, u_copy_args);
  cudaMemcpy(pU->data_d, u, data->numCells * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, u_copy_args);
}

void Poisson::copy_rhs(double *rhs) {
  op_arg rhs_copy_args[] = {
    op_arg_dat(pRHS, -1, OP_ID, 15, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 1, rhs_copy_args);
  cudaMemcpy(rhs, pRHS->data_d, data->numCells * 15 * sizeof(double), cudaMemcpyDeviceToDevice);
  op_mpi_set_dirtybit_cuda(1, rhs_copy_args);
}

void Poisson::create_vec(Vec *v, int size) {
  VecCreateSeqCUDA(PETSC_COMM_SELF, size * data->numCells, v);
}

void Poisson::destroy_vec(Vec *v) {
  VecDestroy(v);
}

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

PetscErrorCode matAMult(Mat A, Vec x, Vec y) {
  Poisson *poisson;
  MatShellGetContext(A, &poisson);
  const double *x_ptr;
  double *y_ptr;
  VecCUDAGetArrayRead(x, &x_ptr);
  VecCUDAGetArray(y, &y_ptr);

  // poisson->rhs(x_ptr, y_ptr);

  VecCUDARestoreArrayRead(x, &x_ptr);
  VecCUDARestoreArray(y, &y_ptr);

  return 0;
}
