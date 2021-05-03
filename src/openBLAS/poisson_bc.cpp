#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_bc(const int numCells, const double *temp0, const double *temp1,
                                double *gx, double *gy) {
  double *mat = (double *)malloc(21 * 15 * sizeof(double));
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, 21, 15, 1.0, constants->invMass, 15, constants->gInterp, 15, 0.0, mat, 15);

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 21, 1.0, mat, 15, temp0, 21, 0.0, gx, 15);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 21, 1.0, mat, 15, temp1, 21, 0.0, gy, 15);
  free(mat);
}

void poisson_bc_blas(INSData *data, Poisson_MF *poisson) {
  // Make sure OP2 data is in the right place
  op_arg poisson_args[] = {
    op_arg_dat(poisson->uFluxX, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->uFluxY, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->gradx, -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(poisson->grady, -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(data->cells, 4, poisson_args);

  openblas_poisson_bc(data->numCells, (double *)poisson->uFluxX->data,
                      (double *)poisson->uFluxY->data, (double *)poisson->gradx->data,
                      (double *)poisson->grady->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(4, poisson_args);
}
