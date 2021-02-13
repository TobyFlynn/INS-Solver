#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_rhs1(const int numCells, const double *fluxXu,
                                  const double *fluxYu, double *qx, double *qy) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, -1.0, LIFT, 15, fluxXu, 15, 1.0, qx, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, -1.0, LIFT, 15, fluxYu, 15, 1.0, qy, 15);
}

void poisson_rhs_blas1(INSData *nsData) {
  // Make sure OP2 data is in the right place
  op_arg poisson_rhs1_args[] = {
    op_arg_dat(nsData->pFluxXu, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->pFluxYu, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->pDuDx, -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(nsData->pDuDy, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges(nsData->cells, 4, poisson_rhs1_args);

  openblas_poisson_rhs1(nsData->numCells, (double *)nsData->pFluxXu->data,
                   (double *)nsData->pFluxYu->data, (double *)nsData->pDuDx->data,
                   (double *)nsData->pDuDy->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(4, poisson_rhs1_args);
}
