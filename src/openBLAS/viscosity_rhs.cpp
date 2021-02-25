#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_viscosity_rhs(const int numCells, const double *qtt0,
                                   const double *qtt1, double *vis0, double *vis1) {
  // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, visMat, 15, qtt0, 15, 0.0, temp0, 15);
  // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, visMat, 15, qtt1, 15, 0.0, temp1, 15);

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, MASS, 15, qtt0, 15, 0.0, vis0, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, MASS, 15, qtt1, 15, 0.0, vis1, 15);
}

void viscosity_rhs_blas(INSData *nsData) {
  // Make sure OP2 data is in the right place
  op_arg viscosity_rhs_args[] = {
    op_arg_dat(nsData->QTT[0], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->QTT[1], -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->visRHS[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->visRHS[1], -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 4, viscosity_rhs_args);

  openblas_viscosity_rhs(nsData->numCells, (double *)nsData->QTT[0]->data,
                      (double *)nsData->QTT[1]->data, (double *)nsData->visRHS[0]->data,
                      (double *)nsData->visRHS[1]->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(4, viscosity_rhs_args);
}
