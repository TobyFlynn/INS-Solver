#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_init_gauss_grad(cublasHandle_t handle, const int numCells,
                        const double *x_d, const double *y_d, double *gxr_d,
                        double *gxs_d, double *gyr_d, double *gys_d) {
  double *g0Dr_d;
  cudaMalloc((void**)&g0Dr_d, 7 * 15 * sizeof(double));
  cudaMemcpy(g0Dr_d, gF0Dr, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  double *g0Ds_d;
  cudaMalloc((void**)&g0Ds_d, 7 * 15 * sizeof(double));
  cudaMemcpy(g0Ds_d, gF0Ds, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  double *g1Dr_d;
  cudaMalloc((void**)&g1Dr_d, 7 * 15 * sizeof(double));
  cudaMemcpy(g1Dr_d, gF1Dr, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  double *g1Ds_d;
  cudaMalloc((void**)&g1Ds_d, 7 * 15 * sizeof(double));
  cudaMemcpy(g1Ds_d, gF1Ds, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  double *g2Dr_d;
  cudaMalloc((void**)&g2Dr_d, 7 * 15 * sizeof(double));
  cudaMemcpy(g2Dr_d, gF2Dr, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);
  double *g2Ds_d;
  cudaMalloc((void**)&g2Ds_d, 7 * 15 * sizeof(double));
  cudaMemcpy(g2Ds_d, gF2Ds, 7 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  // Calc Grad Matrices
  double alpha = 1.0;
  double beta = 0.0;
  for(int c = 0; c < numCells; c++) {
    const double *x = x_d + c * 15;
    const double *y = y_d + c * 15;
    double *gxr = gxr_d + c * 21;
    double *gxs = gxs_d + c * 21;
    double *gyr = gyr_d + c * 21;
    double *gys = gys_d + c * 21;

    // Face 0
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g0Dr_d, 15, x, 1, &beta, gxr, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g0Ds_d, 15, x, 1, &beta, gxs, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g0Dr_d, 15, y, 1, &beta, gyr, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g0Ds_d, 15, y, 1, &beta, gys, 1);

    // Face 1
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g1Dr_d, 15, x, 1, &beta, gxr + 7, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g1Ds_d, 15, x, 1, &beta, gxs + 7, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g1Dr_d, 15, y, 1, &beta, gyr + 7, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g1Ds_d, 15, y, 1, &beta, gys + 7, 1);

    // Face 2
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g2Dr_d, 15, x, 1, &beta, gxr + 14, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g2Ds_d, 15, x, 1, &beta, gxs + 14, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g2Dr_d, 15, y, 1, &beta, gyr + 14, 1);
    cublasDgemv(handle, CUBLAS_OP_T, 15, 7, &alpha, g2Ds_d, 15, y, 1, &beta, gys + 14, 1);
  }

  cudaFree(g0Dr_d);
  cudaFree(g0Ds_d);
  cudaFree(g1Dr_d);
  cudaFree(g1Ds_d);
  cudaFree(g2Dr_d);
  cudaFree(g2Ds_d);
}

void init_gauss_grad_blas(INSData *nsData, GaussData *gaussData) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg init_grad_args[] = {
    op_arg_dat(nsData->x, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->y, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(gaussData->rx, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->sx, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->ry, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->sy, -1, OP_ID, 21, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 6, init_grad_args);

  cublas_init_gauss_grad(handle, nsData->numCells, (double *)nsData->x->data_d,
                   (double *)nsData->y->data_d, (double *)gaussData->rx->data_d,
                   (double *)gaussData->sx->data_d, (double *)gaussData->ry->data_d,
                   (double *)gaussData->sy->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(6, init_grad_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}
