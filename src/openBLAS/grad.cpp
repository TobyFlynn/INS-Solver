#include "cblas.h"

#include "op_seq.h"
#include "../blas_calls.h"

inline void openblas_grad(const int numCells, const double *u, double *div0,
                          double *div1) {
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, Dr, 15, u, 15, 0.0, div0, 15);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 15, numCells, 15, 1.0, Ds, 15, u, 15, 0.0, div1, 15);
}

void grad_blas(INSData *nsData, op_dat u) {
  // Make sure OP2 data is in the right place
  op_arg grad_args[] = {
    op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
    op_arg_dat(nsData->div[0], -1, OP_ID, 15, "double", OP_WRITE),
    op_arg_dat(nsData->div[1], -1, OP_ID, 15, "double", OP_WRITE)
  };
  op_mpi_halo_exchanges(nsData->cells, 3, grad_args);

  openblas_grad(nsData->numCells, (double *)u->data,
               (double *)nsData->div[0]->data, (double *)nsData->div[1]->data);

  // Set correct dirty bits for OP2
  op_mpi_set_dirtybit(3, grad_args);
}
