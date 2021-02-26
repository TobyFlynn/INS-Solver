#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_rhs2(const int numCells, const double *fluxQ,
                                  double *divQ, double *rhs) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, LIFT, 15, fluxQ, 15, -1.0, divQ, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, MASS, 15, divQ, 15, 0.0, rhs, 15);
}

void poisson_rhs_blas2(INSData *nsData, Poisson *pData) {
  // Make sure OP2 data is in the right place
  op_arg poisson_rhs2_args[] = {
    op_arg_dat(pData->pFluxQ, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(pData->pDivQ, -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(pData->pRHS, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 3, poisson_rhs2_args);

  openblas_poisson_rhs2(nsData->numCells, (double *)pData->pFluxQ->data,
                   (double *)pData->pDivQ->data, (double *)pData->pRHS->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(3, poisson_rhs2_args);
}
