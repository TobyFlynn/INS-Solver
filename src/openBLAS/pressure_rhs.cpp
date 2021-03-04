#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_pressure_rhs(const int numCells, const double *div,
                                  const double *dPdN, double *rhs) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, -1.0, MASS, 15, div, 15, 0.0, rhs, 15);
  // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, -1.0, LIFT, 15, dPdN, 15, 1.0, rhs, 15);
}

void pressure_rhs_blas(INSData *nsData, int ind) {
  // Make sure OP2 data is in the right place
  op_arg pressure_rhs_args[] = {
    op_arg_dat(nsData->divVelT, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->dPdN[(ind + 1) % 2], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->pRHS, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 3, pressure_rhs_args);

  openblas_pressure_rhs(nsData->numCells, (double *)nsData->divVelT->data,
                      (double *)nsData->dPdN[(ind + 1) % 2]->data,
                      (double *)nsData->pRHS->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(3, pressure_rhs_args);
}
