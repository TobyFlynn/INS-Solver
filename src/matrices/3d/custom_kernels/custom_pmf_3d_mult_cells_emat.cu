#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

#include "dg_compiler_defs.h"

template<int dg_np, int dg_npf>
__device__ void _pmf_3d_mult_cells_emat_gpu(const int ind, const int *order, const DG_FP *eMat,
                              const DG_FP *in0, const DG_FP *in1, const DG_FP *in2, const DG_FP *in3,
                              DG_FP *out0, DG_FP *out1, DG_FP *out2, DG_FP *out3) {
  const int p = *order;
  const DG_FP *emat_mat = &eMat[(p - 1) * DG_NUM_FACES * DG_NPF * DG_NP];

  out0[ind] = 0.0;
  out1[ind] = 0.0;
  out2[ind] = 0.0;
  out3[ind] = 0.0;
  for(int j = 0; j < DG_NUM_FACES * dg_npf; j++) {
    int mat_ind = DG_MAT_IND(ind, j, dg_np, DG_NUM_FACES * dg_npf);
    out0[ind] += emat_mat[mat_ind] * in0[j];
    out1[ind] += emat_mat[mat_ind] * in1[j];
    out2[ind] += emat_mat[mat_ind] * in2[j];
    out3[ind] += emat_mat[mat_ind] * in3[j];
  }
}

// CUDA kernel function
template<int NUM_CELLS>
__global__ void _op_cuda_pmf_3d_mult_cells_emat(
  const int *__restrict arg0,
  const DG_FP *arg1,
  const DG_FP *__restrict arg2,
  const DG_FP *__restrict arg3,
  const DG_FP *__restrict arg4,
  const DG_FP *__restrict arg5,
  DG_FP *arg6,
  DG_FP *arg7,
  DG_FP *arg8,
  DG_FP *arg9,
  int   set_size ) {

  __shared__ DG_FP emat_shared[DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP];
  __shared__ DG_FP in0_shared[NUM_CELLS * DG_NUM_FACES * DG_NPF];
  __shared__ DG_FP in1_shared[NUM_CELLS * DG_NUM_FACES * DG_NPF];
  __shared__ DG_FP in2_shared[NUM_CELLS * DG_NUM_FACES * DG_NPF];
  __shared__ DG_FP in3_shared[NUM_CELLS * DG_NUM_FACES * DG_NPF];

  for(int i = threadIdx.x; i < DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP; i += blockDim.x) {
    emat_shared[i] = arg1[i];
  }

  __syncthreads();

  for (int n = threadIdx.x + blockIdx.x * blockDim.x; n < set_size * DG_NP; n += blockDim.x * gridDim.x){
    const int node_id = n % DG_NP;
    const int cell_id = n / DG_NP;
    const int local_cell_id = (n / DG_NP) - ((n - threadIdx.x) / DG_NP);

    // If entire thread is in set
    if(n - threadIdx.x + blockDim.x < set_size * DG_NP) {
      __syncthreads();
      const int start_ind = ((n - threadIdx.x) / DG_NP) * DG_NUM_FACES * DG_NPF;
      const int num_elem  = ((n - threadIdx.x + blockDim.x) / DG_NP) - ((n - threadIdx.x) / DG_NP) + 1;
      for(int i = threadIdx.x; i < num_elem * DG_NUM_FACES * DG_NPF; i += blockDim.x) {
        in0_shared[i] = arg2[start_ind + i];
        in1_shared[i] = arg3[start_ind + i];
        in2_shared[i] = arg4[start_ind + i];
        in3_shared[i] = arg5[start_ind + i];
      }
      __syncthreads();

      switch(*(arg0 + cell_id * 1)) {
        case 1:
          _pmf_3d_mult_cells_emat_gpu<4, 3>(node_id, arg0 + cell_id * 1,
                                 emat_shared,
                                 in0_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 in1_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 in2_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 in3_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 arg6 + cell_id * DG_NP,
                                 arg7 + cell_id * DG_NP,
                                 arg8 + cell_id * DG_NP,
                                 arg9 + cell_id * DG_NP);
          break;
        case 2:
          _pmf_3d_mult_cells_emat_gpu<10, 6>(node_id, arg0 + cell_id * 1,
                                 emat_shared,
                                 in0_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 in1_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 in2_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 in3_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 arg6 + cell_id * DG_NP,
                                 arg7 + cell_id * DG_NP,
                                 arg8 + cell_id * DG_NP,
                                 arg9 + cell_id * DG_NP);
          break;
        case 3:
          _pmf_3d_mult_cells_emat_gpu<20, 10>(node_id, arg0 + cell_id * 1,
                                 emat_shared,
                                 in0_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 in1_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 in2_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 in3_shared + local_cell_id * DG_NUM_FACES * DG_NPF,
                                 arg6 + cell_id * DG_NP,
                                 arg7 + cell_id * DG_NP,
                                 arg8 + cell_id * DG_NP,
                                 arg9 + cell_id * DG_NP);
          break;
      }
    } else {
      switch(*(arg0 + cell_id * 1)) {
        case 1:
          _pmf_3d_mult_cells_emat_gpu<4, 3>(node_id, arg0 + cell_id * 1,
                                 emat_shared,
                                 arg2 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg3 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg4 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg5 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg6 + cell_id * DG_NP,
                                 arg7 + cell_id * DG_NP,
                                 arg8 + cell_id * DG_NP,
                                 arg9 + cell_id * DG_NP);
          break;
        case 2:
          _pmf_3d_mult_cells_emat_gpu<10, 6>(node_id, arg0 + cell_id * 1,
                                 emat_shared,
                                 arg2 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg3 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg4 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg5 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg6 + cell_id * DG_NP,
                                 arg7 + cell_id * DG_NP,
                                 arg8 + cell_id * DG_NP,
                                 arg9 + cell_id * DG_NP);
          break;
        case 3:
          _pmf_3d_mult_cells_emat_gpu<20, 10>(node_id, arg0 + cell_id * 1,
                                 emat_shared,
                                 arg2 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg3 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg4 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg5 + cell_id * DG_NUM_FACES * DG_NPF,
                                 arg6 + cell_id * DG_NP,
                                 arg7 + cell_id * DG_NP,
                                 arg8 + cell_id * DG_NP,
                                 arg9 + cell_id * DG_NP);
          break;
      }
    }
  }
}


//host stub function
void custom_kernel_pmf_3d_mult_cells_emat(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  DG_FP*arg1h = (DG_FP *)arg1.data;
  int nargs = 10;
  op_arg args[10];

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


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pmf_3d_mult_cells_emat");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP * sizeof(DG_FP));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg1.data   = OP_consts_h + consts_bytes;
    arg1.data_d = OP_consts_d + consts_bytes;
    memcpy(arg1.data, arg1h, DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP * sizeof(DG_FP));
    // for ( int d=0; d<2400; d++ ){
    //   ((DG_FP *)arg1.data)[d] = arg1h[d];
    // }
    consts_bytes += ROUND_UP(DG_ORDER * DG_NUM_FACES * DG_NPF * DG_NP * sizeof(DG_FP));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    const int nthread = 256;
    const int nblocks = 200 < (set->size * DG_NP) / nthread + 1 ? 200 : (set->size * DG_NP) / nthread + 1;
    const int num_cells = (nthread / DG_NP) + 2;

    _op_cuda_pmf_3d_mult_cells_emat<num_cells><<<nblocks,nthread>>>(
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
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
