#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_set_rhs(const int numCells, double *liftx, double *lifty, double *tau) {
  double *tempx = (double *)malloc(numCells * 15 * sizeof(double));
  double *tempy = (double *)malloc(numCells * 15 * sizeof(double));
  double *tempt = (double *)malloc(numCells * 15 * sizeof(double));

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, LIFT, 15, liftx, 15, 0.0, tempx, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, LIFT, 15, lifty, 15, 0.0, tempy, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, LIFT, 15, tau, 15, 0.0, tempt, 15);

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, MASS, 15, tempx, 15, 0.0, liftx, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, MASS, 15, tempy, 15, 0.0, lifty, 15);

  // memcpy(liftx, tempx, numCells * 15 * sizeof(double));
  // memcpy(lifty, tempy, numCells * 15 * sizeof(double));
  memcpy(tau, tempt, numCells * 15 * sizeof(double));

  free(tempx);
  free(tempy);
  free(tempt);
}

void poisson_set_rhs_blas(INSData *nsData, Poisson *pData) {
  // Make sure OP2 data is in the right place
  op_arg set_rhs_args[] = {
    op_arg_dat(pData->pExRHS[0], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(pData->pExRHS[1], -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(pData->bcTau, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges(nsData->cells, 3, set_rhs_args);

  openblas_poisson_set_rhs(nsData->numCells, (double *)pData->pExRHS[0]->data,
                           (double *)pData->pExRHS[1]->data, (double *)pData->bcTau->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(3, set_rhs_args);
}
