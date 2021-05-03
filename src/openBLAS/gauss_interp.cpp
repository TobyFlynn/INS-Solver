#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_gauss_interp(const int numCells, const double *in, double *out) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 21, numCells, 15, 1.0, constants->gInterp, 15, in, 15, 0.0, out, 21);
}

void gauss_interp_blas(INSData *data, op_dat input, op_dat output) {
  // Make sure OP2 data is in the right place
  op_arg interp_args[] = {
    op_arg_dat(input, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(output, -1, OP_ID, 21, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 2, interp_args);

  openblas_gauss_interp(data->numCells, (double *)input->data, (double *)output->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(2, interp_args);
}
