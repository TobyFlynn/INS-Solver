#include "op_seq.h"

#include <memory>
#include <iostream>
#include <stdexcept>

#include "dg_utils.h"

DG_FP *getOP2PtrDevice(op_dat dat, op_access acc) {
  throw std::runtime_error("\ngetOP2PtrDevice not implemented for CPU\n");
}

void releaseOP2PtrDevice(op_dat dat, op_access acc, const DG_FP *ptr) {
  throw std::runtime_error("\releaseOP2PtrDevice not implemented for CPU\n");
}

DG_FP *getOP2PtrHost(op_dat dat, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_halo_exchanges(dat->set, 1, args);
  op_mpi_wait_all(1, args);
  return (DG_FP *)dat->data;
}

void releaseOP2PtrHost(op_dat dat, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };

  op_mpi_set_dirtybit(1, args);

  ptr = nullptr;
}

DG_FP *getOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc) {
  throw std::runtime_error("\ngetOP2PtrDevice not implemented for CPU\n");
}

void releaseOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc, const DG_FP *ptr) {
  throw std::runtime_error("\releaseOP2PtrDevice not implemented for CPU\n");
}

DG_FP *getOP2PtrHostMap(op_dat dat, op_map map, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, DG_FP_STR, acc),
    op_arg_dat(dat, 1, map, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_halo_exchanges(map->from, 2, args);
  op_mpi_wait_all(2, args);
  return (DG_FP *)dat->data;
}

void releaseOP2PtrHostMap(op_dat dat, op_map map, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, DG_FP_STR, acc),
    op_arg_dat(dat, 1, map, dat->dim, DG_FP_STR, acc)
  };

  op_mpi_set_dirtybit(2, args);

  ptr = nullptr;
}
