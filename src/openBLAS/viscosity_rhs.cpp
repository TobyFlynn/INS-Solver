#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_viscosity_rhs(const int numCells, double *qtt0, double *qtt1) {
  double *temp0 = (double *)malloc(15 * numCells * sizeof(double));
  double *temp1 = (double *)malloc(15 * numCells * sizeof(double));

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, visMat, 15, qtt0, 15, 0.0, temp0, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, visMat, 15, qtt1, 15, 0.0, temp1, 15);

  memcpy(qtt0, temp0, 15 * numCells * sizeof(double));
  memcpy(qtt1, temp1, 15 * numCells * sizeof(double));

  free(temp0);
  free(temp1);
}

void viscosity_rhs_blas(INSData *nsData) {
  // Make sure OP2 data is in the right place
  op_arg viscosity_rhs_args[] = {
    op_arg_dat(nsData->QTT[0], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(nsData->QTT[1], -1, OP_ID, 15, "double", OP_RW),
  };
  op_mpi_halo_exchanges(nsData->cells, 2, viscosity_rhs_args);

  openblas_viscosity_rhs(nsData->numCells, (double *)nsData->QTT[0]->data,
                      (double *)nsData->QTT[1]->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(2, viscosity_rhs_args);
}
