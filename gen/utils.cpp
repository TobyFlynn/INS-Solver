#include "op_seq.h"

#include <memory>

double *getOP2Array(op_dat dat) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, 16, "double", OP_READ)
  };
  op_mpi_halo_exchanges(dat->set, 1, args);
  double *res = (double *)malloc(dat->set->size * 16 * sizeof(double));
  memcpy(res, dat->data, dat->set->size * 16 * sizeof(double));
  op_mpi_set_dirtybit(1, args);
  return res;
}
