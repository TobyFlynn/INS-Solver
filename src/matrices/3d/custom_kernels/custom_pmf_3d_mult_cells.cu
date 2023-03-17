#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

#include "dg_compiler_defs.h"

template<int dg_np>
__device__ void shared_pmf_3d_mult_cells_gpu(const int ind, const int *p, const DG_FP *dr,
                            const DG_FP *ds, const DG_FP *dt,
                            const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                            const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                            const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                            const DG_FP *in_x, const DG_FP *in_y,
                            const DG_FP *in_z, DG_FP *out) {
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];

  for(int n = 0; n < dg_np; n++) {
    int mat_ind = DG_MAT_IND(n, ind, dg_np, dg_np);
    out[ind] += dr_mat[mat_ind] * (rx[0] * in_x[n] + ry[0] * in_y[n] + rz[0] * in_z[n]);
    out[ind] += ds_mat[mat_ind] * (sx[0] * in_x[n] + sy[0] * in_y[n] + sz[0] * in_z[n]);
    out[ind] += dt_mat[mat_ind] * (tx[0] * in_x[n] + ty[0] * in_y[n] + tz[0] * in_z[n]);
  }
}

template<int dg_np>
__device__ void _pmf_3d_mult_cells_gpu(const int ind, const int *p, const DG_FP *dr,
                            const DG_FP *ds, const DG_FP *dt,
                            const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                            const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                            const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                            const DG_FP *in_x, const DG_FP *in_y,
                            const DG_FP *in_z, DG_FP *out) {
  const DG_FP *dr_mat = &dr[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *ds_mat = &ds[(*p - 1) * DG_NP * DG_NP];
  const DG_FP *dt_mat = &dt[(*p - 1) * DG_NP * DG_NP];

  for(int n = 0; n < dg_np; n++) {
    int mat_ind = DG_MAT_IND(n, ind, dg_np, dg_np);
    out[ind] += dr_mat[mat_ind] * (rx[0] * in_x[n] + ry[0] * in_y[n] + rz[0] * in_z[n]);
    out[ind] += ds_mat[mat_ind] * (sx[0] * in_x[n] + sy[0] * in_y[n] + sz[0] * in_z[n]);
    out[ind] += dt_mat[mat_ind] * (tx[0] * in_x[n] + ty[0] * in_y[n] + tz[0] * in_z[n]);
  }
}

// CUDA kernel function
template<int NUM_CELLS>
__global__ void _op_cuda_pmf_3d_mult_cells(
  const int *__restrict arg0,
  const DG_FP *arg1,
  const DG_FP *arg2,
  const DG_FP *arg3,
  const DG_FP *__restrict arg4,
  const DG_FP *__restrict arg5,
  const DG_FP *__restrict arg6,
  const DG_FP *__restrict arg7,
  const DG_FP *__restrict arg8,
  const DG_FP *__restrict arg9,
  const DG_FP *__restrict arg10,
  const DG_FP *__restrict arg11,
  const DG_FP *__restrict arg12,
  const DG_FP *__restrict arg16,
  const DG_FP *__restrict arg17,
  const DG_FP *__restrict arg18,
  DG_FP *arg19,
  int   set_size ) {

  // Load matrices into shared memory
  __shared__ DG_FP dr_shared[DG_ORDER * DG_NP * DG_NP];
  __shared__ DG_FP ds_shared[DG_ORDER * DG_NP * DG_NP];
  __shared__ DG_FP dt_shared[DG_ORDER * DG_NP * DG_NP];
  __shared__ DG_FP ux_shared[NUM_CELLS * DG_NP];
  __shared__ DG_FP uy_shared[NUM_CELLS * DG_NP];
  __shared__ DG_FP uz_shared[NUM_CELLS * DG_NP];

  for(int i = threadIdx.x; i < DG_ORDER * DG_NP * DG_NP; i += blockDim.x) {
    dr_shared[i] = arg1[i];
    ds_shared[i] = arg2[i];
    dt_shared[i] = arg3[i];
  }

  __syncthreads();

  for (int n = threadIdx.x + blockIdx.x * blockDim.x; n < set_size * DG_NP; n += blockDim.x * gridDim.x){
    const int node_id = n % DG_NP;
    const int cell_id = n / DG_NP;
    const int local_cell_id = (n / DG_NP) - ((n - threadIdx.x) / DG_NP);

    // If entire thread is in set
    if(n - threadIdx.x + blockDim.x < set_size * DG_NP) {
      __syncthreads();
      const int start_ind = ((n - threadIdx.x) / DG_NP) * DG_NP;
      const int num_elem  = ((n - threadIdx.x + blockDim.x) / DG_NP) - ((n - threadIdx.x) / DG_NP) + 1;
      for(int i = threadIdx.x; i < num_elem * DG_NP; i += blockDim.x) {
        ux_shared[i] = arg16[start_ind + i];
        uy_shared[i] = arg17[start_ind + i];
        uz_shared[i] = arg18[start_ind + i];
      }

      __syncthreads();

      switch(*(arg0 + cell_id * 1)) {
        case 1:
          shared_pmf_3d_mult_cells_gpu<4>(node_id, arg0 + cell_id * 1,
                            dr_shared,
                            ds_shared,
                            dt_shared,
                            arg4 + cell_id * 1,
                            arg5 + cell_id * 1,
                            arg6 + cell_id * 1,
                            arg7 + cell_id * 1,
                            arg8 + cell_id * 1,
                            arg9 + cell_id * 1,
                            arg10 + cell_id * 1,
                            arg11 + cell_id * 1,
                            arg12 + cell_id * 1,
                            ux_shared + local_cell_id * DG_NP,
                            uy_shared + local_cell_id * DG_NP,
                            uz_shared + local_cell_id * DG_NP,
                            arg19 + cell_id * DG_NP);
          break;
        case 2:
          shared_pmf_3d_mult_cells_gpu<10>(node_id, arg0 + cell_id * 1,
                            dr_shared,
                            ds_shared,
                            dt_shared,
                            arg4 + cell_id * 1,
                            arg5 + cell_id * 1,
                            arg6 + cell_id * 1,
                            arg7 + cell_id * 1,
                            arg8 + cell_id * 1,
                            arg9 + cell_id * 1,
                            arg10 + cell_id * 1,
                            arg11 + cell_id * 1,
                            arg12 + cell_id * 1,
                            ux_shared + local_cell_id * DG_NP,
                            uy_shared + local_cell_id * DG_NP,
                            uz_shared + local_cell_id * DG_NP,
                            arg19 + cell_id * DG_NP);
          break;
        case 3:
          shared_pmf_3d_mult_cells_gpu<20>(node_id, arg0 + cell_id * 1,
                            dr_shared,
                            ds_shared,
                            dt_shared,
                            arg4 + cell_id * 1,
                            arg5 + cell_id * 1,
                            arg6 + cell_id * 1,
                            arg7 + cell_id * 1,
                            arg8 + cell_id * 1,
                            arg9 + cell_id * 1,
                            arg10 + cell_id * 1,
                            arg11 + cell_id * 1,
                            arg12 + cell_id * 1,
                            ux_shared + local_cell_id * DG_NP,
                            uy_shared + local_cell_id * DG_NP,
                            uz_shared + local_cell_id * DG_NP,
                            arg19 + cell_id * DG_NP);
          break;
      }
    } else {
      switch(*(arg0 + cell_id * 1)) {
        case 1:
          _pmf_3d_mult_cells_gpu<4>(node_id, arg0 + cell_id * 1,
                            dr_shared,
                            ds_shared,
                            dt_shared,
                            arg4 + cell_id * 1,
                            arg5 + cell_id * 1,
                            arg6 + cell_id * 1,
                            arg7 + cell_id * 1,
                            arg8 + cell_id * 1,
                            arg9 + cell_id * 1,
                            arg10 + cell_id * 1,
                            arg11 + cell_id * 1,
                            arg12 + cell_id * 1,
                            arg16 + cell_id * DG_NP,
                            arg17 + cell_id * DG_NP,
                            arg18 + cell_id * DG_NP,
                            arg19 + cell_id * DG_NP);
          break;
        case 2:
          _pmf_3d_mult_cells_gpu<10>(node_id, arg0 + cell_id * 1,
                            dr_shared,
                            ds_shared,
                            dt_shared,
                            arg4 + cell_id * 1,
                            arg5 + cell_id * 1,
                            arg6 + cell_id * 1,
                            arg7 + cell_id * 1,
                            arg8 + cell_id * 1,
                            arg9 + cell_id * 1,
                            arg10 + cell_id * 1,
                            arg11 + cell_id * 1,
                            arg12 + cell_id * 1,
                            arg16 + cell_id * DG_NP,
                            arg17 + cell_id * DG_NP,
                            arg18 + cell_id * DG_NP,
                            arg19 + cell_id * DG_NP);
          break;
        case 3:
          _pmf_3d_mult_cells_gpu<20>(node_id, arg0 + cell_id * 1,
                            dr_shared,
                            ds_shared,
                            dt_shared,
                            arg4 + cell_id * 1,
                            arg5 + cell_id * 1,
                            arg6 + cell_id * 1,
                            arg7 + cell_id * 1,
                            arg8 + cell_id * 1,
                            arg9 + cell_id * 1,
                            arg10 + cell_id * 1,
                            arg11 + cell_id * 1,
                            arg12 + cell_id * 1,
                            arg16 + cell_id * DG_NP,
                            arg17 + cell_id * DG_NP,
                            arg18 + cell_id * DG_NP,
                            arg19 + cell_id * DG_NP);
          break;
      }
    }
  }
/*
  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    switch(*(arg0+n*1)) {
      case 1:
        _pmf_3d_mult_cells_gpu<4>(arg0+n*1,
                          dr_shared,
                          ds_shared,
                          dt_shared,
                          arg4+n*1,
                          arg5+n*1,
                          arg6+n*1,
                          arg7+n*1,
                          arg8+n*1,
                          arg9+n*1,
                          arg10+n*1,
                          arg11+n*1,
                          arg12+n*1,
                          arg13+n*DG_NP,
                          arg14+n*DG_NP,
                          arg15+n*DG_NP,
                          arg16+n*DG_NP,
                          arg17+n*DG_NP,
                          arg18+n*DG_NP,
                          arg19+n*DG_NP);
        break;
      case 2:
        _pmf_3d_mult_cells_gpu<10>(arg0+n*1,
                          dr_shared,
                          ds_shared,
                          dt_shared,
                          arg4+n*1,
                          arg5+n*1,
                          arg6+n*1,
                          arg7+n*1,
                          arg8+n*1,
                          arg9+n*1,
                          arg10+n*1,
                          arg11+n*1,
                          arg12+n*1,
                          arg13+n*DG_NP,
                          arg14+n*DG_NP,
                          arg15+n*DG_NP,
                          arg16+n*DG_NP,
                          arg17+n*DG_NP,
                          arg18+n*DG_NP,
                          arg19+n*DG_NP);
        break;
      case 3:
        _pmf_3d_mult_cells_gpu<20>(arg0+n*1,
                          dr_shared,
                          ds_shared,
                          dt_shared,
                          arg4+n*1,
                          arg5+n*1,
                          arg6+n*1,
                          arg7+n*1,
                          arg8+n*1,
                          arg9+n*1,
                          arg10+n*1,
                          arg11+n*1,
                          arg12+n*1,
                          arg13+n*DG_NP,
                          arg14+n*DG_NP,
                          arg15+n*DG_NP,
                          arg16+n*DG_NP,
                          arg17+n*DG_NP,
                          arg18+n*DG_NP,
                          arg19+n*DG_NP);
        break;
    }
  }
  */
}


//host stub function
void custom_kernel_pmf_3d_mult_cells(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16){

  DG_FP*arg1h = (DG_FP *)arg1.data;
  DG_FP*arg2h = (DG_FP *)arg2.data;
  DG_FP*arg3h = (DG_FP *)arg3.data;
  int nargs = 17;
  op_arg args[17];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  args[10] = arg10;
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;
  args[16] = arg16;

  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pmf_3d_mult_cells");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(DG_ORDER * DG_NP * DG_NP * sizeof(DG_FP));
    consts_bytes += ROUND_UP(DG_ORDER * DG_NP * DG_NP * sizeof(DG_FP));
    consts_bytes += ROUND_UP(DG_ORDER * DG_NP * DG_NP * sizeof(DG_FP));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg1.data   = OP_consts_h + consts_bytes;
    arg1.data_d = OP_consts_d + consts_bytes;
    memcpy(arg1.data, arg1h, DG_ORDER * DG_NP * DG_NP * sizeof(DG_FP));
    consts_bytes += ROUND_UP(DG_ORDER * DG_NP * DG_NP * sizeof(DG_FP));
    arg2.data   = OP_consts_h + consts_bytes;
    arg2.data_d = OP_consts_d + consts_bytes;
    memcpy(arg2.data, arg2h, DG_ORDER * DG_NP * DG_NP * sizeof(DG_FP));
    consts_bytes += ROUND_UP(DG_ORDER * DG_NP * DG_NP * sizeof(DG_FP));
    arg3.data   = OP_consts_h + consts_bytes;
    arg3.data_d = OP_consts_d + consts_bytes;
    memcpy(arg3.data, arg3h, DG_ORDER * DG_NP * DG_NP * sizeof(DG_FP));
    consts_bytes += ROUND_UP(DG_ORDER * DG_NP * DG_NP * sizeof(DG_FP));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    const int nthread = 256;
    const int nblocks = 200 < (set->size * DG_NP) / nthread + 1 ? 200 : (set->size * DG_NP) / nthread + 1;
    const int num_cells = (nthread / DG_NP) + 2;

    _op_cuda_pmf_3d_mult_cells<num_cells><<<nblocks,nthread>>>(
      (int *) arg0.data_d,
      (DG_FP *) arg1.data_d,
      (DG_FP *) arg2.data_d,
      (DG_FP *) arg3.data_d,
      (DG_FP *) arg4.data_d,
      (DG_FP *) arg5.data_d,
      (DG_FP *) arg6.data_d,
      (DG_FP *) arg7.data_d,
      (DG_FP *) arg8.data_d,
      (DG_FP *) arg9.data_d,
      (DG_FP *) arg10.data_d,
      (DG_FP *) arg11.data_d,
      (DG_FP *) arg12.data_d,
      (DG_FP *) arg13.data_d,
      (DG_FP *) arg14.data_d,
      (DG_FP *) arg15.data_d,
      (DG_FP *) arg16.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
