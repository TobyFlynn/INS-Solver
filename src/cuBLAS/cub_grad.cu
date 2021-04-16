#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_cub_grad(cublasHandle_t handle, const int numCells,
                            const double *in_d, double *out_d) {
  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 46, numCells, 15, &alpha, constants->cubV_d, 15, in_d, 15, &beta, out_d, 46);
}

void cub_grad_blas(INSData *data, CubatureData *cubatureData, op_dat u) {
  // Make sure OP2 data is in the right place
  op_arg v_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(cubatureData->op_temps[0], -1, OP_ID, 46, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(data->cells, 2, v_args);

  cublas_cub_grad(constants->handle, data->numCells, (double *)u->data_d,
                  (double *)cubatureData->op_temps[0]->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(2, v_args);
}
