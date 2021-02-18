#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_viscosity_rhs(cublasHandle_t handle, const int numCells,
                                 double *qtt0_d, double *qtt1_d) {
  double *visMat_d;
  cudaMalloc((void**)&visMat_d, 15 * 15 * sizeof(double));
  cudaMemcpy(visMat_d, visMat, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *temp0_d;
  cudaMalloc((void**)&temp0_d, 15 * numCells * sizeof(double));
  double *temp1_d;
  cudaMalloc((void**)&temp1_d, 15 * numCells * sizeof(double));

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, visMat_d, 15, qtt0_d, 15, &beta, temp0_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, visMat_d, 15, qtt1_d, 15, &beta, temp1_d, 15);

  cudaMemcpy(qtt0_d, temp0_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(qtt1_d, temp1_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

  cudaFree(temp0_d);
  cudaFree(temp1_d);
  cudaFree(visMat_d);
}

void viscosity_rhs_blas(INSData *nsData) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg viscosity_rhs_args[] = {
    op_arg_dat(nsData->QTT[0], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(nsData->QTT[1], -1, OP_ID, 15, "double", OP_RW),
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 2, viscosity_rhs_args);

  cublas_viscosity_rhs(handle, nsData->numCells, (double *)nsData->QTT[0]->data_d,
                      (double *)nsData->QTT[1]->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(2, viscosity_rhs_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
