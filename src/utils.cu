#include "op_seq.h"

#include <memory>

#include "dg_utils.h"

double *getOP2PtrDevice(op_dat dat, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, "double", acc)
  };
  op_mpi_halo_exchanges_cuda(dat->set, 1, args);
  op_mpi_wait_all(1, args);
  return (double *) dat->data_d;
}

void releaseOP2PtrDevice(op_dat dat, op_access acc, const double *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, "double", acc)
  };
  op_mpi_set_dirtybit_cuda(1, args);

  ptr = nullptr;
}

double *getOP2PtrHost(op_dat dat, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, "double", acc)
  };
  op_mpi_halo_exchanges_cuda(dat->set, 1, args);
  op_mpi_wait_all(1, args);
  double *res = (double *)malloc(dat->set->size * dat->dim * sizeof(double));
  cudaMemcpy(res, dat->data_d, dat->set->size * dat->dim * sizeof(double), cudaMemcpyDeviceToHost);
  return res;
}

void releaseOP2PtrHost(op_dat dat, op_access acc, const double *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, "double", acc)
  };

  if(acc != OP_READ) {
    cudaMemcpy(dat->data_d, ptr, dat->set->size * dat->dim * sizeof(double), cudaMemcpyHostToDevice);
  }

  op_mpi_set_dirtybit_cuda(1, args);

  free((void *)ptr);
  ptr = nullptr;
}

double *getOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, "double", acc),
    op_arg_dat(dat, 1, map, dat->dim, "double", acc)
  };
  op_mpi_halo_exchanges_cuda(map->from, 2, args);
  op_mpi_wait_all(2, args);
  
  return (double *) dat->data_d;
}

void releaseOP2PtrDeviceMap(op_dat dat, op_map map, op_access acc, const double *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, "double", acc),
    op_arg_dat(dat, 1, map, dat->dim, "double", acc)
  };
  op_mpi_set_dirtybit_cuda(2, args);

  ptr = nullptr;
}

double *getOP2PtrHostMap(op_dat dat, op_map map, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, "double", acc),
    op_arg_dat(dat, 1, map, dat->dim, "double", acc)
  };
  op_mpi_halo_exchanges_cuda(map->from, 2, args);
  op_mpi_wait_all(2, args);

  // TODO check that this shouldn't be more than dat->set->size
  double *res = (double *)malloc(dat->set->size * dat->dim * sizeof(double));
  cudaMemcpy(res, dat->data_d, dat->set->size * dat->dim * sizeof(double), cudaMemcpyDeviceToHost);
  return res;
}

void releaseOP2PtrHostMap(op_dat dat, op_map map, op_access acc, const double *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, 0, map, dat->dim, "double", acc),
    op_arg_dat(dat, 1, map, dat->dim, "double", acc)
  };

  if(acc != OP_READ) {
    cudaMemcpy(dat->data_d, ptr, dat->set->size * dat->dim * sizeof(double), cudaMemcpyHostToDevice);
  }

  op_mpi_set_dirtybit(2, args);
  free((void *)ptr);
  ptr = nullptr;
}

double *getOP2Array(op_dat dat) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, DG_SUB_CELLS, "double", OP_READ)
  };
  op_mpi_halo_exchanges_cuda(dat->set, 1, args);
  double *res = (double *)malloc(dat->set->size * DG_SUB_CELLS * sizeof(double));
  cudaMemcpy(res, dat->data_d, dat->set->size * DG_SUB_CELLS * sizeof(double), cudaMemcpyDeviceToHost);
  op_mpi_set_dirtybit_cuda(1, args);
  return res;
}
