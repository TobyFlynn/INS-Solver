#include "cublas_v2.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void cublas_div(cublasHandle_t handle, const int numCells, const double *u_d,
                       const double *v_d, double *div0_d, double *div1_d,
                       double *div2_d, double *div3_d) {
  // CUBLAS_OP_T because cublas is column major but constants are stored row major
  double alpha = 1.0;
  double beta = 0.0;
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, constants->Dr_d, 15, u_d, 15, &beta, div0_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, constants->Ds_d, 15, u_d, 15, &beta, div1_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, constants->Dr_d, 15, v_d, 15, &beta, div2_d, 15);
  cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, 15, numCells, 15, &alpha, constants->Ds_d, 15, v_d, 15, &beta, div3_d, 15);
}

void div_blas(INSData *nsData, op_dat u, op_dat v) {
  // Make sure OP2 data is in the right place
  op_arg div_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(v, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->div[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->div[1], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->div[2], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->div[3], -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges_cuda(nsData->cells, 6, div_args);

  cublas_div(constants->handle, nsData->numCells, (double *)u->data_d,
                   (double *)v->data_d, (double *)nsData->div[0]->data_d,
                   (double *)nsData->div[1]->data_d, (double *)nsData->div[2]->data_d,
                   (double *)nsData->div[3]->data_d);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit_cuda(6, div_args);
}
