#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_mass(const int numCells, const double *pU,
                                  double *pRHS, double factor) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, factor, MASS, 15, pU, 15, 1.0, pRHS, 15);
}

void poisson_mass_blas(INSData *nsData, Poisson *pData, double factor) {
  // Make sure OP2 data is in the right place
  op_arg poisson_mass_args[] = {
    op_arg_dat(pData->pU, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(pData->pRHS, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges(nsData->cells, 2, poisson_mass_args);

  openblas_poisson_mass(nsData->numCells, (double *)pData->pU->data,
                   (double *)pData->pRHS->data, factor);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(2, poisson_mass_args);
}
