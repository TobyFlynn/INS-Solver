#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_poisson_test_rhs(cublasHandle_t handle, const int numCells,
                                    double *rhs_d) {
  double *temp_d;
  cudaMalloc((void**)&temp_d, 15 * numCells * sizeof(double));

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, constants->mass_d, 15, rhs_d, 15, &beta, temp_d, 15);

  cudaMemcpy(rhs_d, temp_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

  cudaFree(temp_d);
}

void poisson_test_rhs_blas(INSData *nsData, op_dat rhs) {
  // Make sure OP2 data is in the right place
  op_arg poisson_test_rhs_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 1, poisson_test_rhs_args);

  cublas_poisson_test_rhs(constants->handle, nsData->numCells, (double *)rhs->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(1, poisson_test_rhs_args);
}
