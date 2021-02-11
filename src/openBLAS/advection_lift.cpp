#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_advection_lift(const int numCells, const double *flux0,
                               const double *flux1, double *N0, double *N1) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, LIFT, 15, flux0, 15, 1.0, N0, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, LIFT, 15, flux1, 15, 1.0, N1, 15);
}

void advection_lift_blas(INSData *nsData) {
  // Make sure OP2 data is in the right place
  op_arg advec_lift_args[] = {
    op_arg_dat(nsData->flux[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->flux[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->N[0], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(nsData->N[1], -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges(nsData->cells, 4, advec_lift_args);

  openblas_advection_lift(nsData->numCells, (double *)nsData->flux[0]->data,
                     (double *)nsData->flux[1]->data, (double *)nsData->N[0]->data,
                     (double *)nsData->N[1]->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(4, advec_lift_args);
}
