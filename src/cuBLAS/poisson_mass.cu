#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_poisson_mass(cublasHandle_t handle, const int numCells,
                        const double *pU_d, double *pRHS_d, double factor) {
  double *MASS_d;
  cudaMalloc((void**)&MASS_d, 15 * 15 * sizeof(double));
  cudaMemcpy(MASS_d, MASS, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = factor;
  double beta = 1.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, MASS_d, 15, pU_d, 15, &beta, pRHS_d, 15);

  cudaFree(MASS_d);
}

void poisson_mass_blas(INSData *nsData, Poisson *pData, double factor) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg poisson_mass_args[] = {
    op_arg_dat(pData->pU, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(pData->pRHS, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 2, poisson_mass_args);

  cublas_poisson_mass(handle, nsData->numCells, (double *)pData->pU->data_d,
                   (double *)pData->pRHS->data_d, factor);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(2, poisson_mass_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
