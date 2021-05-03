#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_gauss_coords(const int numCells, const double *x, const double *y,
                                  double *gx, double *gy) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 21, numCells, 15, 1.0, constants->gInterp, 15, x, 15, 0.0, gx, 21);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 21, numCells, 15, 1.0, constants->gInterp, 15, y, 15, 0.0, gy, 21);
}

void init_gauss_coords_blas(INSData *nsData, GaussData *gaussData) {
  // Make sure OP2 data is in the right place
  op_arg gauss_args[] = {
    op_arg_dat(nsData->x, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->y, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(gaussData->x, -1, OP_ID, 21, "double", OP_WRITE),
    op_arg_dat(gaussData->y, -1, OP_ID, 21, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 4, gauss_args);

  openblas_gauss_coords(nsData->numCells, (double *)nsData->x->data,
                        (double *)nsData->y->data, (double *)gaussData->x->data,
                        (double *)gaussData->y->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(4, gauss_args);
}
