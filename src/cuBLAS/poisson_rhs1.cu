#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_poisson_rhs1(cublasHandle_t handle, const int numCells,
                                const double *fluxX_d, const double *fluxY_d,
                                double *qx_d, double *qy_d) {
  double *interp_d;
  cudaMalloc((void**)&interp_d, 21 * 15 * sizeof(double));
  cudaMemcpy(interp_d, gInterp, 21 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *invMass_d;
  cudaMalloc((void**)&invMass_d, 15 * 15 * sizeof(double));
  cudaMemcpy(invMass_d, invMass, 15 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *temp_d;
  cudaMalloc((void**)&temp_d, 15 * numCells * sizeof(double));

  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = -1.0;
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 21, &alpha, interp_d, 15, fluxX_d, 21, &beta, qx_d, 15);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 21, &alpha, interp_d, 15, fluxY_d, 21, &beta, qy_d, 15);

  double beta2 = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, invMass_d, 15, qx_d, 15, &beta2, temp_d, 15);
  cudaMemcpy(qx_d, temp_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, invMass_d, 15, qy_d, 15, &beta2, temp_d, 15);
  cudaMemcpy(qy_d, temp_d, 15 * numCells * sizeof(double), cudaMemcpyDeviceToDevice);

  cudaFree(interp_d);
  cudaFree(invMass_d);
  cudaFree(temp_d);
}

void poisson_rhs_blas1(INSData *data, Poisson_MF *poisson) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg poisson_args[] = {
    op_arg_dat(poisson->uFluxX, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->uFluxY, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->qx, -1, OP_ID, 21, "double", OP_RW),
    op_arg_dat(poisson->qy, -1, OP_ID, 21, "double", OP_RW)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 4, poisson_args);

  cublas_poisson_rhs1(handle, data->numCells, (double *)poisson->uFluxX->data_d,
                      (double *)poisson->uFluxY->data_d, (double *)poisson->qx->data_d,
                      (double *)poisson->qy->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(4, poisson_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
