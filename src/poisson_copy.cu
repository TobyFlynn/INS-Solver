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
