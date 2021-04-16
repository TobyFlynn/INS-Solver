#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_init_cubature(const int numCells, const double *x,
                                const double *y, double *cxr, double *cxs,
                                double *cyr, double *cys) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 46, numCells, 15, 1.0, constants->cubDr, 15, x, 15, 0.0, cxr, 46);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 46, numCells, 15, 1.0, constants->cubDs, 15, x, 15, 0.0, cxs, 46);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 46, numCells, 15, 1.0, constants->cubDr, 15, y, 15, 0.0, cyr, 46);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 46, numCells, 15, 1.0, constants->cubDs, 15, y, 15, 0.0, cys, 46);
}

void init_cubature_blas(INSData *nsData, CubatureData *cubData) {
  // Make sure OP2 data is in the right place
  op_arg init_cubature_args[] = {
    op_arg_dat(nsData->x, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->y, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(cubData->rx, -1, OP_ID, 46, "double", OP_WRITE),
    op_arg_dat(cubData->sx, -1, OP_ID, 46, "double", OP_WRITE),
    op_arg_dat(cubData->ry, -1, OP_ID, 46, "double", OP_WRITE),
    op_arg_dat(cubData->sy, -1, OP_ID, 46, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 6, init_cubature_args);

  openblas_init_cubature(nsData->numCells, (double *)nsData->x->data,
                     (double *)nsData->y->data, (double *)cubData->rx->data,
                     (double *)cubData->sx->data, (double *)cubData->ry->data,
                     (double *)cubData->sy->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(6, init_cubature_args);
}
