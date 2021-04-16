#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_poisson_test_rhs(const int numCells, double *rhs) {
  double *temp = (double *)malloc(15 * numCells * sizeof(double));

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, constants->MASS, 15, rhs, 15, 0.0, temp, 15);

  memcpy(rhs, temp, 15 * numCells * sizeof(double));

  free(temp);
}

void poisson_test_rhs_blas(INSData *nsData, op_dat rhs) {
  // Make sure OP2 data is in the right place
  op_arg poisson_test_rhs_args[] = {
    op_arg_dat(rhs, -1, OP_ID, 15, "double", OP_RW)
  };
  op_mpi_halo_exchanges(nsData->cells, 1, poisson_test_rhs_args);

  openblas_poisson_test_rhs(nsData->numCells, (double *)rhs->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(1, poisson_test_rhs_args);
}
