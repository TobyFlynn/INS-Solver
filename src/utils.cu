#include "op_seq.h"

#ifdef INS_MPI
#include "op_lib_mpi.h"
#endif

#include <memory>

#include "dg_utils.h"

DG_FP *getOP2PtrDevice(op_dat dat, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_halo_exchanges_grouped(dat->set, 1, args, 2);
  op_mpi_wait_all_grouped(1, args, 2);
  return (DG_FP *) dat->data_d;
}

void releaseOP2PtrDevice(op_dat dat, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_set_dirtybit_cuda(1, args);

  ptr = nullptr;
}

DG_FP *getOP2PtrHost(op_dat dat, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_halo_exchanges_grouped(dat->set, 1, args, 2);
  op_mpi_wait_all_grouped(1, args, 2);
  const int size = getSetSizeFromOpArg(&args[0]);
  DG_FP *res = (DG_FP *)malloc(size * dat->dim * sizeof(DG_FP));
  cudaMemcpy(res, dat->data_d, size * dat->dim * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  return res;
}

void releaseOP2PtrHost(op_dat dat, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };

  if(acc != OP_READ) {
    const int size = getSetSizeFromOpArg(&args[0]);
    cudaMemcpy(dat->data_d, ptr, size * dat->dim * sizeof(DG_FP), cudaMemcpyHostToDevice);
  }

  op_mpi_set_dirtybit_cuda(1, args);

  free((void *)ptr);
  ptr = nullptr;
}

DG_FP *getOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, DG_FP_STR, acc),
    op_arg_dat(dat, 1, map, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_halo_exchanges_grouped(map->from, 2, args, 2);
  op_mpi_wait_all_grouped(2, args, 2);

  return (DG_FP *) dat->data_d;
}

void releaseOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, DG_FP_STR, acc),
    op_arg_dat(dat, 1, map, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_set_dirtybit_cuda(2, args);

  ptr = nullptr;
}

DG_FP *getOP2PtrHostMap(op_dat dat, op_map map, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, DG_FP_STR, acc),
    op_arg_dat(dat, 1, map, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_halo_exchanges_grouped(map->from, 2, args, 2);
  op_mpi_wait_all_grouped(2, args, 2);

  const int size = getSetSizeFromOpArg(&args[0]);
  DG_FP *res = (DG_FP *)malloc(dat->dim * size * sizeof(DG_FP));
  cudaMemcpy(res, dat->data_d, size * dat->dim * sizeof(DG_FP), cudaMemcpyDeviceToHost);

  return res;
}

void releaseOP2PtrHostMap(op_dat dat, op_map map, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, DG_FP_STR, acc),
    op_arg_dat(dat, 1, map, dat->dim, DG_FP_STR, acc)
  };

  if(acc != OP_READ) {
    const int size = getSetSizeFromOpArg(&args[0]);
    cudaMemcpy(dat->data_d, ptr, size * dat->dim * sizeof(DG_FP), cudaMemcpyHostToDevice);
  }

  op_mpi_set_dirtybit_cuda(2, args);

  free((void *)ptr);
  ptr = nullptr;
}
