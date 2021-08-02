#include "op_seq.h"

#include <memory>

double *getOP2Array(op_dat dat) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, 9, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(dat->set, 1, args);
  double *res = (double *)malloc(dat->set->size * 9 * sizeof(double));
  cudaMemcpy(res, dat->data_d, dat->set->size * 9 * sizeof(double), cudaMemcpyDeviceToHost);
  op_mpi_set_dirtybit_cuda(1, args);
  return res;
}
