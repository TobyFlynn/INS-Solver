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

#if DG_DIM == 3
#include "dg_global_constants/dg_mat_constants_dev_ptrs_3d.h"
#include "op_cuda_rt_support.h"

__constant__ DG_FP *dg_r_kernel;
__constant__ DG_FP *dg_s_kernel;
__constant__ DG_FP *dg_t_kernel;
__constant__ DG_FP *dg_Dr_kernel;
__constant__ DG_FP *dg_Ds_kernel;
__constant__ DG_FP *dg_Dt_kernel;
__constant__ DG_FP *dg_Drw_kernel;
__constant__ DG_FP *dg_Dsw_kernel;
__constant__ DG_FP *dg_Dtw_kernel;
__constant__ DG_FP *dg_Mass_kernel;
__constant__ DG_FP *dg_InvMass_kernel;
__constant__ DG_FP *dg_InvV_kernel;
__constant__ DG_FP *dg_Lift_kernel;
__constant__ DG_FP *dg_MM_F0_kernel;
__constant__ DG_FP *dg_MM_F1_kernel;
__constant__ DG_FP *dg_MM_F2_kernel;
__constant__ DG_FP *dg_MM_F3_kernel;
__constant__ DG_FP *dg_Emat_kernel;
__constant__ DG_FP *dg_Interp_kernel;

void transfer_kernel_ptrs() {
  cutilSafeCall(cudaMemcpyToSymbol(dg_r_kernel, &dg_r_d, sizeof(dg_r_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_s_kernel, &dg_s_d, sizeof(dg_s_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_t_kernel, &dg_t_d, sizeof(dg_t_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Dr_kernel, &dg_Dr_d, sizeof(dg_Dr_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Ds_kernel, &dg_Ds_d, sizeof(dg_Ds_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Dt_kernel, &dg_Dt_d, sizeof(dg_Dt_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Drw_kernel, &dg_Drw_d, sizeof(dg_Drw_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Dsw_kernel, &dg_Dsw_d, sizeof(dg_Dsw_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Dtw_kernel, &dg_Dtw_d, sizeof(dg_Dtw_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Mass_kernel, &dg_Mass_d, sizeof(dg_Mass_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_InvMass_kernel, &dg_InvMass_d, sizeof(dg_InvMass_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Lift_kernel, &dg_Lift_d, sizeof(dg_Lift_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_MM_F0_kernel, &dg_MM_F0_d, sizeof(dg_MM_F0_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_MM_F1_kernel, &dg_MM_F1_d, sizeof(dg_MM_F1_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_MM_F2_kernel, &dg_MM_F2_d, sizeof(dg_MM_F2_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_MM_F3_kernel, &dg_MM_F3_d, sizeof(dg_MM_F3_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Emat_kernel, &dg_Emat_d, sizeof(dg_Emat_d)));
  cutilSafeCall(cudaMemcpyToSymbol(dg_Interp_kernel, &dg_Interp_d, sizeof(dg_Interp_d)));
}
#endif
