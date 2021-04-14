#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_poisson_test_rhs(cublasHandle_t handle, const int numCells,
                                    double *rhs_d) {
  double *MASS_d;
  cudaMalloc((void**)&MASS_d, 15 * 15 * sizeof(double));
  cudaMemcpy(MASS_d, MASS, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *temp_d;
  cudaMalloc((void**)&temp_d, 15 * numCells * sizeof(double));

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, MASS_d, 15, rhs_d, 15, &beta, temp_d, 15);

  cudaMemcpy(rhs_d, temp_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

  cudaFree(temp_d);
  cudaFree(MASS_d);
}

void poisson_test_rhs_blas(INSData *nsData, op_dat rhs) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg poisson_test_rhs_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 1, poisson_test_rhs_args);

  cublas_poisson_test_rhs(handle, nsData->numCells, (double *)rhs->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(1, poisson_test_rhs_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
