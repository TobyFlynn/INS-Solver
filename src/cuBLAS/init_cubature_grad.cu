#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_init_cubature(cublasHandle_t handle, const int numCells,
                        const double *x_d, const double *y_d, double *cxr_d,
                        double *cxs_d, double *cyr_d, double *cys_d) {
  double *cubVDr_d;
  cudaMalloc((void**)&cubVDr_d, 46 * 15 * sizeof(double));
  cudaMemcpy(cubVDr_d, cubVDr, 46 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double *cubVDs_d;
  cudaMalloc((void**)&cubVDs_d, 46 * 15 * sizeof(double));
  cudaMemcpy(cubVDs_d, cubVDs, 46 * 15 * sizeof(double), cudaMemcpyHostToDevice);

  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 46, numCells, 15, &alpha, cubVDr_d, 15, x_d, 15, &beta, cxr_d, 46);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 46, numCells, 15, &alpha, cubVDs_d, 15, x_d, 15, &beta, cxs_d, 46);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 46, numCells, 15, &alpha, cubVDr_d, 15, y_d, 15, &beta, cyr_d, 46);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 46, numCells, 15, &alpha, cubVDs_d, 15, y_d, 15, &beta, cys_d, 46);

  cudaFree(cubVDr_d);
  cudaFree(cubVDs_d);
}

void init_cubature_grad_blas(INSData *nsData, CubatureData *cubData) {
  // Initialise cuBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);
  cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST);
  // Make sure OP2 data is in the right place
  op_arg init_cubature_args[] = {
    op_arg_dat(nsData->x, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->y, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(cubData->rx, -1, OP_ID, 46, "double", OP_WRITE),
    op_arg_dat(cubData->sx, -1, OP_ID, 46, "double", OP_WRITE),
    op_arg_dat(cubData->ry, -1, OP_ID, 46, "double", OP_WRITE),
    op_arg_dat(cubData->sy, -1, OP_ID, 46, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 6, init_cubature_args);

  cublas_init_cubature(handle, nsData->numCells, (double *)nsData->x->data_d,
                   (double *)nsData->y->data_d, (double *)cubData->rx->data_d,
                   (double *)cubData->sx->data_d, (double *)cubData->ry->data_d,
                   (double *)cubData->sy->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(6, init_cubature_args);
  // Free resources used by cuBLAS
  cublasDestroy(handle);
}