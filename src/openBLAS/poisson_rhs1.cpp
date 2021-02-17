#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_rhs1(const int numCells, const double *fluxXu,
                                  const double *fluxYu, double *qx, double *qy) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, -1.0, LIFT, 15, fluxXu, 15, 1.0, qx, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, -1.0, LIFT, 15, fluxYu, 15, 1.0, qy, 15);
}

void poisson_rhs_blas1(INSData *nsData, Poisson *pData) {
  // Make sure OP2 data is in the right place
  op_arg poisson_rhs1_args[] = {
    op_arg_dat(pData->pFluxXu, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(pData->pFluxYu, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(pData->pDuDx, -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(pData->pDuDy, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges(nsData->cells, 4, poisson_rhs1_args);

  openblas_poisson_rhs1(nsData->numCells, (double *)pData->pFluxXu->data,
                   (double *)pData->pFluxYu->data, (double *)pData->pDuDx->data,
                   (double *)pData->pDuDy->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(4, poisson_rhs1_args);
}
