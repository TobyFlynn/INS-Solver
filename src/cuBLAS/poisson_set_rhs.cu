#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_poisson_set_rhs(cublasHandle_t handle, const int numCells,
                                   double *liftx_d, double *lifty_d, double *tau_d) {
  double *LIFT_d;
  cudaMalloc((void**)&LIFT_d, 15 * 15 * sizeof(double));
  cudaMemcpy(LIFT_d, LIFT, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *MASS_d;
  cudaMalloc((void**)&MASS_d, 15 * 15 * sizeof(double));
  cudaMemcpy(MASS_d, MASS, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *tempx_d;
  cudaMalloc((void**)&tempx_d, 15 * numCells * sizeof(double));
  double *tempy_d;
  cudaMalloc((void**)&tempy_d, 15 * numCells * sizeof(double));
  double *tempt_d;
  cudaMalloc((void**)&tempy_d, 15 * numCells * sizeof(double));

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, LIFT_d, 15, liftx_d, 15, &beta, tempx_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, LIFT_d, 15, lifty_d, 15, &beta, tempy_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, LIFT_d, 15, tau_d, 15, &beta, tempt_d, 15);

  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, MASS_d, 15, tempx_d, 15, &beta, liftx_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, MASS_d, 15, tempy_d, 15, &beta, lifty_d, 15);

  // cudaMemcpy(liftx_d, tempx_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  // cudaMemcpy(lifty_d, tempy_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(tau_d, tempt_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

  cudaFree(LIFT_d);
  cudaFree(MASS_d);
  cudaFree(tempx_d);
  cudaFree(tempy_d);
  cudaFree(tempt_d);
}

void poisson_set_rhs_blas(INSData *nsData, Poisson *pData) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg set_rhs_args[] = {
    op_arg_dat(pData->pExRHS[0], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(pData->pExRHS[1], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(pData->bcTau, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 3, set_rhs_args);

  cublas_poisson_set_rhs(handle, nsData->numCells, (double *)pData->pExRHS[0]->data_d,
                         (double *)pData->pExRHS[1]->data_d, (double *)pData->bcTau->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(3, set_rhs_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
