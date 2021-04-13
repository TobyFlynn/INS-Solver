#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_rhs1(const int numCells, const double *fluxX,
                                  const double *fluxY, double *qx, double *qy) {
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 21, 1.0, gInterp, 15, fluxX, 21, -1.0, qx, 15);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 21, 1.0, gInterp, 15, fluxY, 21, -1.0, qy, 15);

  double *temp = (double *)malloc(15 * numCells * sizeof(double));
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, invMass, 15, qx, 15, 0.0, temp, 15);
  memcpy(qx, temp, 15 * numCells * sizeof(double));

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, invMass, 15, qy, 15, 0.0, temp, 15);
  memcpy(qy, temp, 15 * numCells * sizeof(double));

  free(temp);
}

void poisson_rhs_blas1(INSData *data, Poisson_MF *poisson) {
  // Make sure OP2 data is in the right place
  op_arg poisson_args[] = {
    op_arg_dat(poisson->uFluxX, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->uFluxY, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->qx, -1, OP_ID, 21, "double", OP_RW),
    op_arg_dat(poisson->qy, -1, OP_ID, 21, "double", OP_RW)
  };
  op_mpi_halo_exchanges(data->cells, 4, poisson_args);

  openblas_poisson_rhs1(data->numCells, (double *)poisson->uFluxX->data,
                        (double *)poisson->uFluxY->data, (double *)poisson->qx->data,
                        (double *)poisson->qy->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(4, poisson_args);
}
