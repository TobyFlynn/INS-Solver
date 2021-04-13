#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_rhs2(const int numCells, const double *flux, double *rhs) {
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 15, numCells, 21, -1.0, gInterp, 15, flux, 21, 1.0, rhs, 15);
}

void poisson_rhs_blas2(INSData *data, Poisson_MF *poisson) {
  // Make sure OP2 data is in the right place
  op_arg poisson_args[] = {
    op_arg_dat(poisson->qFlux, -1, OP_ID, 21, "double", OP_READ),
    op_arg_dat(poisson->rhs, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges(data->cells, 2, poisson_args);

  openblas_poisson_rhs2(data->numCells, (double *)poisson->qFlux->data,
                        (double *)poisson->rhs->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(2, poisson_args);
}
