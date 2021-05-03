#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_cub_grad_w(const int numCells, const double *in, double *out) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 46, numCells, 15, 1.0, constants->cubV, 15, in, 15, 0.0, out, 46);
}

void cub_grad_w_blas(INSData *data, CubatureData *cubatureData, op_dat u) {
  // Make sure OP2 data is in the right place
  op_arg v_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(cubatureData->op_temps[0], -1, OP_ID, 46, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 2, v_args);

  openblas_cub_grad_w(data->numCells, (double *)u->data, (double *)cubatureData->op_temps[0]->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(2, v_args);
}
