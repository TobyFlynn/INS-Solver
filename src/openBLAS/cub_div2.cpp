#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_cub_div2(const int numCells, const double *temp0, double *res) {
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 46, 1.0, constants->cubV, 15, temp0, 46, 0.0, res, 15);
}

void cub_div_blas2(INSData *data, CubatureData *cubatureData, op_dat res) {
  // Make sure OP2 data is in the right place
  op_arg div_args[] = {
    op_arg_dat(cubatureData->op_temps[0], -1, OP_ID, 46, "double", OP_READ),
    op_arg_dat(res, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 2, div_args);

  openblas_cub_div2(data->numCells, (double *)cubatureData->op_temps[0]->data,
                    (double *)res->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(2, div_args);
}
