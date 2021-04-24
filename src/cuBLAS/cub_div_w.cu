#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_cub_div_w(cublasHandle_t handle, const int numCells,
                           const double *in0_d, const double *in1_d,
                           double *out0_d, double *out1_d) {
  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  cudaStream_t stream0, stream1;
  cudaStreamCreate(&stream0);
  cudaStreamCreate(&stream1);

  double alpha = 1.0;
  double beta = 0.0;
  cublasSetStream(handle, stream0);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 46, numCells, 15, &alpha, constants->cubV_d, 15, in0_d, 15, &beta, out0_d, 46);
  cublasSetStream(handle, stream1);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 46, numCells, 15, &alpha, constants->cubV_d, 15, in1_d, 15, &beta, out1_d, 46);

  cudaStreamDestroy(stream0);
  cudaStreamDestroy(stream1);
  cublasSetStream(handle, NULL);
}

void cub_div_w_blas(INSData *data, CubatureData *cubatureData, op_dat u, op_dat v) {
  // Make sure OP2 data is in the right place
  op_arg v_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(v, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(cubatureData->op_temps[0], -1, OP_ID, 46, "double", OP_WRITE),
    op_arg_dat(cubatureData->op_temps[1], -1, OP_ID, 46, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 4, v_args);

  cublas_cub_div_w(constants->handle, data->numCells, (double *)u->data_d, (double *)v->data_d,
                  (double *)cubatureData->op_temps[0]->data_d, (double *)cubatureData->op_temps[1]->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(4, v_args);
}