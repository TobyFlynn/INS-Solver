#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_cub_grad2_w(cublasHandle_t handle, const int numCells,
                            const double *temp0, const double *temp1,
                            const double *temp2, const double *temp3,
                            double *ux_d, double *uy_d) {
  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  cudaStream_t stream0, stream1;
  cudaStreamCreate(&stream0);
  cudaStreamCreate(&stream1);

  double alpha = 1.0;
  double beta = 0.0;
  double beta2 = 1.0;
  cublasSetStream(handle, stream0);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 46, &alpha, constants->cubDr_d, 15, temp0, 46, &beta, ux_d, 15);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 46, &alpha, constants->cubDs_d, 15, temp1, 46, &beta2, ux_d, 15);
  cublasSetStream(handle, stream1);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 46, &alpha, constants->cubDr_d, 15, temp2, 46, &beta, uy_d, 15);
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 15, numCells, 46, &alpha, constants->cubDs_d, 15, temp3, 46, &beta2, uy_d, 15);

  cudaStreamDestroy(stream0);
  cudaStreamDestroy(stream1);
  cublasSetStream(handle, NULL);
}

void cub_grad_w_blas2(INSData *data, CubatureData *cubatureData, op_dat ux, op_dat uy) {
  // Make sure OP2 data is in the right place
  op_arg grad_args[] = {
    op_arg_dat(cubatureData->op_temps[0], -1, OP_ID, 46, "double", OP_READ),
    op_arg_dat(cubatureData->op_temps[1], -1, OP_ID, 46, "double", OP_READ),
    op_arg_dat(cubatureData->op_temps[2], -1, OP_ID, 46, "double", OP_READ),
    op_arg_dat(cubatureData->op_temps[3], -1, OP_ID, 46, "double", OP_READ),
    op_arg_dat(ux, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(uy, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 6, grad_args);

  cublas_cub_grad2_w(constants->handle, data->numCells, (double *)cubatureData->op_temps[0]->data_d,
                  (double *)cubatureData->op_temps[1]->data_d, (double *)cubatureData->op_temps[2]->data_d,
                  (double *)cubatureData->op_temps[3]->data_d, (double *)ux->data_d, (double *)uy->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(6, grad_args);
}
