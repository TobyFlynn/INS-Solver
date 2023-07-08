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
  op_mpi_halo_exchanges_grouped(dat->set, 1, args, 2, 0);
  op_mpi_wait_all_grouped(1, args, 2, 0);
  return (DG_FP *) dat->data_d;
}

void releaseOP2PtrDevice(op_dat dat, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_set_dirtybit_cuda(1, args);

  ptr = nullptr;
}

__global__ void soa_to_aos_utils(const int set_size, const int stride,
                                 const int dim, const DG_FP *in, DG_FP *out) {
  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= set_size * dim) return;
  const int node = tid / set_size;
  const int cell = tid % set_size;
  const int in_ind = cell + node * stride;
  const int out_ind = cell * dim + node;
  out[out_ind] = in[in_ind];
}

DG_FP *getOP2PtrHost(op_dat dat, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_halo_exchanges_grouped(dat->set, 1, args, 2, 0);
  op_mpi_wait_all_grouped(1, args, 2, 0);

  const int size = getSetSizeFromOpArg(&args[0]);
  DG_FP *res = (DG_FP *)malloc(size * dat->dim * sizeof(DG_FP));
  #ifdef DG_OP2_SOA
  DG_FP *res_d;
  cudaMalloc(&res_d, size * dat->dim * sizeof(DG_FP));
  const int nthread = 512;
  const int nblocks = (dat->set->size * dat->dim - 1) / nthread + 1;
  soa_to_aos_utils<<<nblocks,nthread>>>(dat->set->size, size, dat->dim, (DG_FP *)dat->data_d, res_d);
  cudaMemcpy(res, res_d, size * dat->dim * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  cudaFree(res_d);
  #else
  cudaMemcpy(res, dat->data_d, size * dat->dim * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  #endif
  return res;
}

__global__ void aos_to_soa_utils(const int set_size, const int stride,
                                 const int dim, const DG_FP *in, DG_FP *out) {
  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= set_size * dim) return;
  const int in_ind = tid;
  const int node = tid % dim;
  const int cell = tid / dim;
  const int out_ind = cell + node * stride;
  out[out_ind] = in[in_ind];
}

void releaseOP2PtrHost(op_dat dat, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };

  if(acc != OP_READ) {
    const int size = getSetSizeFromOpArg(&args[0]);
    #ifdef DG_OP2_SOA
    DG_FP *ptr_d;
    cudaMalloc(&ptr_d, size * dat->dim * sizeof(DG_FP));
    cudaMemcpy(ptr_d, ptr, size * dat->dim * sizeof(DG_FP), cudaMemcpyHostToDevice);
    const int nthread = 512;
    const int nblocks = (dat->set->size * dat->dim - 1) / nthread + 1;
    aos_to_soa_utils<<<nblocks,nthread>>>(dat->set->size, size, dat->dim, ptr_d, (DG_FP *)dat->data_d);
    cudaFree(ptr_d);
    #else
    cudaMemcpy(dat->data_d, ptr, size * dat->dim * sizeof(DG_FP), cudaMemcpyHostToDevice);
    #endif
  }

  op_mpi_set_dirtybit_cuda(1, args);

  free((void *)ptr);
  ptr = nullptr;
}

DG_FP *getOP2PtrDeviceHE(op_dat dat, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc),
  };
  op_mpi_halo_exchanges_grouped(dat->set, 1, args, 2, 1);
  op_mpi_wait_all_grouped(1, args, 2, 1);

  return (DG_FP *) dat->data_d;
}

void releaseOP2PtrDeviceHE(op_dat dat, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc),
  };
  op_mpi_set_dirtybit_cuda(1, args);

  ptr = nullptr;
}

DG_FP *getOP2PtrHostHE(op_dat dat, op_access acc) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };
  op_mpi_halo_exchanges_grouped(dat->set, 1, args, 2, 1);
  op_mpi_wait_all_grouped(1, args, 2, 1);

  const int size = getSetSizeFromOpArg(&args[0]);
  const int set_size = dat->set->size + dat->set->exec_size + dat->set->nonexec_size;
  DG_FP *res = (DG_FP *)malloc(size * dat->dim * sizeof(DG_FP));
  #ifdef DG_OP2_SOA
  DG_FP *res_d;
  cudaMalloc(&res_d, size * dat->dim * sizeof(DG_FP));
  const int nthread = 512;
  const int nblocks = (set_size * dat->dim - 1) / nthread + 1;
  soa_to_aos_utils<<<nblocks,nthread>>>(set_size, size, dat->dim, (DG_FP *)dat->data_d, res_d);
  cudaMemcpy(res, res_d, size * dat->dim * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  cudaFree(res_d);
  #else
  cudaMemcpy(res, dat->data_d, size * dat->dim * sizeof(DG_FP), cudaMemcpyDeviceToHost);
  #endif
  return res;
}

void releaseOP2PtrHostHE(op_dat dat, op_access acc, const DG_FP *ptr) {
  op_arg args[] = {
    op_arg_dat(dat, -1, OP_ID, dat->dim, DG_FP_STR, acc)
  };

  if(acc != OP_READ) {
    const int size = getSetSizeFromOpArg(&args[0]);
    const int set_size = dat->set->size + dat->set->exec_size + dat->set->nonexec_size;
    #ifdef DG_OP2_SOA
    DG_FP *ptr_d;
    cudaMalloc(&ptr_d, size * dat->dim * sizeof(DG_FP));
    cudaMemcpy(ptr_d, ptr, size * dat->dim * sizeof(DG_FP), cudaMemcpyHostToDevice);
    const int nthread = 512;
    const int nblocks = (set_size * dat->dim - 1) / nthread + 1;
    aos_to_soa_utils<<<nblocks,nthread>>>(set_size, size, dat->dim, ptr_d, (DG_FP *)dat->data_d);
    cudaFree(ptr_d);
    #else
    cudaMemcpy(dat->data_d, ptr, size * dat->dim * sizeof(DG_FP), cudaMemcpyHostToDevice);
    #endif
  }

  op_mpi_set_dirtybit_cuda(1, args);

  free((void *)ptr);
  ptr = nullptr;
}
