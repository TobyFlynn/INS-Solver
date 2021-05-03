#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_rhs1(const int numCells, const double *fluxX,
                                  const double *fluxY, double *dudx, double *dudy,
                                  double *qx, double *qy, double *gradx, double *grady) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, constants->invMass, 15, dudx, 15, 0.0, gradx, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, constants->invMass, 15, dudy, 15, 0.0, grady, 15);

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 21, -1.0, constants->gInterp, 15, fluxX, 21, 1.0, dudx, 15);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 21, -1.0, constants->gInterp, 15, fluxY, 21, 1.0, dudy, 15);

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, constants->invMass, 15, dudx, 15, 0.0, qx, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, constants->invMass, 15, dudy, 15, 0.0, qy, 15);
}

void poisson_rhs_blas1(INSData *data, Poisson_MF *poisson) {
  // Make sure OP2 data is in the right place
  op_arg poisson_args[] = {
    op_arg_dat(poisson->uFluxX, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->uFluxY, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->dudx, -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(poisson->dudy, -1, OP_ID, 15, "double", OP_RW),
    op_arg_dat(poisson->qx, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(poisson->qy, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(poisson->gradx, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(poisson->grady, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 8, poisson_args);

  openblas_poisson_rhs1(data->numCells, (double *)poisson->uFluxX->data,
                        (double *)poisson->uFluxY->data, (double *)poisson->dudx->data,
                        (double *)poisson->dudy->data, (double *)poisson->qx->data,
                        (double *)poisson->qy->data, (double *)poisson->gradx->data,
                        (double *)poisson->grady->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(8, poisson_args);
}
