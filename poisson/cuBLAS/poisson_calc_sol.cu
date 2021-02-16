#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_poisson_calc_sol(cublasHandle_t handle, const int numCells,
                                   double *sol_d) {
  double *invM_d;
  cudaMalloc((void**)&invM_d, 15 * 15 * sizeof(double));
  cudaMemcpy(invM_d, invM, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *temp_d;
  cudaMalloc((void**)&temp_d, 15 * numCells * sizeof(double));

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, invM_d, 15, sol_d, 15, &beta, temp_d, 15);

  cudaMemcpy(sol_d, temp_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

  cudaFree(temp_d);
  cudaFree(invM_d);
}

void poisson_calc_sol_blas(INSData *nsData) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg poisson_sol_args[] = {
    op_arg_dat(nsData->sol, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 1, poisson_sol_args);

  cublas_poisson_calc_sol(handle, nsData->numCells, (double *)nsData->sol->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(1, poisson_sol_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
