#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

#include "dg_compiler_defs.h"
#include "dg_global_constants/dg_mat_constants_3d.h"

template<int dg_np, int NS>
__device__ void _fpmf_grad_3d_gpu(const int ind, const DG_FP *dr, const DG_FP *ds,
                             const DG_FP *dt, const DG_FP *u, const DG_FP *fact,
                             const DG_FP *geof, DG_FP *ux, DG_FP *uy, DG_FP *uz) {
  if(!(ind < dg_np)) return;

  DG_FP tmp_r[NS] = {};
  DG_FP tmp_s[NS] = {};
  DG_FP tmp_t[NS] = {};
  #pragma unroll
  for(int n = 0; n < dg_np; n++) {
    const DG_FP u_n = u[n];
    #pragma unroll
    for(int i = 0; i < NS; i++) {
      int mat_ind = DG_MAT_IND(ind + i, n, dg_np, dg_np);
      tmp_r[i] += dr[mat_ind] * u_n;
      tmp_s[i] += ds[mat_ind] * u_n;
      tmp_t[i] += dt[mat_ind] * u_n;
    }
  }

  #pragma unroll
  for(int i = 0; i < NS; i++) {
    const DG_FP fact_ = fact[ind + i];
    ux[ind + i] = fact_ * (geof[RX_IND] * tmp_r[i] + geof[SX_IND] * tmp_s[i] + geof[TX_IND] * tmp_t[i]);
    uy[ind + i] = fact_ * (geof[RY_IND] * tmp_r[i] + geof[SY_IND] * tmp_s[i] + geof[TY_IND] * tmp_t[i]);
    uz[ind + i] = fact_ * (geof[RZ_IND] * tmp_r[i] + geof[SZ_IND] * tmp_s[i] + geof[TZ_IND] * tmp_t[i]);
  }
}

// CUDA kernel function
template<int p, int NUM_CELLS, int NS>
__global__ void _op_cuda_fpmf_grad_3d(
  const int *__restrict arg0,
  const DG_FP *__restrict arg4,
  const DG_FP *__restrict argFactor,
  const DG_FP *__restrict argGeof,
  DG_FP *arg14,
  DG_FP *arg15,
  DG_FP *arg16,
  int   set_size ) {
  const int np = (p + 1) * (p + 2) * (p + 3) / 6;
  // Load matrices into shared memory
  __shared__ DG_FP dr_shared[DG_NP * DG_NP];
  __shared__ DG_FP ds_shared[DG_NP * DG_NP];
  __shared__ DG_FP dt_shared[DG_NP * DG_NP];
  __shared__ DG_FP u_shared[NS * NUM_CELLS * DG_NP];

  const int start_ind_mat = (p - 1) * DG_NP * DG_NP;
  for(int i = threadIdx.x; i < DG_NP * DG_NP; i += blockDim.x) {
    dr_shared[i] = dg_Dr_kernel[start_ind_mat + i];
    ds_shared[i] = dg_Ds_kernel[start_ind_mat + i];
    dt_shared[i] = dg_Dt_kernel[start_ind_mat + i];
  }

  __syncthreads();

  const int n = threadIdx.x * NS + blockIdx.x * blockDim.x * NS;
  if(n < set_size * DG_NP) {
    const int node_id = n % DG_NP;
    const int cell_id = n / DG_NP;
    const int local_cell_id = (n / DG_NP) - ((n - NS * threadIdx.x) / DG_NP);
    // If entire thread is in set
    if(n - NS * threadIdx.x + NS * blockDim.x < set_size * DG_NP) {
      const int start_ind = ((n - NS * threadIdx.x) / DG_NP) * DG_NP;
      const int num_elem  = ((n - NS * threadIdx.x + NS * blockDim.x) / DG_NP) - ((n - NS * threadIdx.x) / DG_NP) + 1;
      for(int i = threadIdx.x; i < num_elem * DG_NP; i += blockDim.x) {
        u_shared[i] = arg4[start_ind + i];
      }
      __syncthreads();

      _fpmf_grad_3d_gpu<np,NS>(node_id,
                dr_shared,
                ds_shared,
                dt_shared,
                u_shared + local_cell_id * DG_NP, //arg4 + cell_id * DG_NP,
                argFactor + cell_id * DG_NP,
                argGeof + cell_id * 10,
                arg14 + cell_id * DG_NP,
                arg15 + cell_id * DG_NP,
                arg16 + cell_id * DG_NP);
    } else {
      _fpmf_grad_3d_gpu<np,NS>(node_id,
                dr_shared,
                ds_shared,
                dt_shared,
                arg4 + cell_id * DG_NP,
                argFactor + cell_id * DG_NP,
                argGeof + cell_id * 10,
                arg14 + cell_id * DG_NP,
                arg15 + cell_id * DG_NP,
                arg16 + cell_id * DG_NP);
    }
  }
}

#include "timing.h"
extern Timing *timer;

//host stub function
void custom_kernel_fpmf_grad_3d(const int order, char const *name, op_set set,
  op_arg arg0,
  op_arg arg4,
  op_arg argFactor,
  op_arg argGeof,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16){

  int nargs = 7;
  op_arg args[7];

  args[0] = arg0;
  args[1] = arg4;
  args[2] = argFactor;
  args[3] = argGeof;
  args[4] = arg14;
  args[5] = arg15;
  args[6] = arg16;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  grad_3d");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {
    //set CUDA execution parameters
    const int nthread = 256;
    // const int nblocks = 200 < (set->size * DG_NP) / nthread + 1 ? 200 : (set->size * DG_NP) / nthread + 1;
    const int nblocks = order == 3 ? (set->size * DG_NP) / (nthread * 2) + 1 : (set->size * DG_NP) / nthread + 1;
    const int num_cells = (nthread / DG_NP) + 2;

    cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_grad_3d<1,num_cells,1>, cudaFuncCachePreferShared));
    cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_grad_3d<2,num_cells,1>, cudaFuncCachePreferShared));
    cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_grad_3d<3,num_cells,2>, cudaFuncCachePreferShared));
    cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_grad_3d<4,num_cells,1>, cudaFuncCachePreferShared));
    cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_grad_3d<5,num_cells,1>, cudaFuncCachePreferShared));

    switch(order) {
      case 1:
        _op_cuda_fpmf_grad_3d<1,num_cells,1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
          (DG_FP *) argGeof.data_d,
          (DG_FP *) arg14.data_d,
          (DG_FP *) arg15.data_d,
          (DG_FP *) arg16.data_d,
          set->size );
        break;
      case 2:
        _op_cuda_fpmf_grad_3d<2,num_cells,1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
          (DG_FP *) argGeof.data_d,
          (DG_FP *) arg14.data_d,
          (DG_FP *) arg15.data_d,
          (DG_FP *) arg16.data_d,
          set->size );
        break;
      case 3:
        timer->startTimer("fpmf_grad 3rd order");
        _op_cuda_fpmf_grad_3d<3,num_cells,2><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
          (DG_FP *) argGeof.data_d,
          (DG_FP *) arg14.data_d,
          (DG_FP *) arg15.data_d,
          (DG_FP *) arg16.data_d,
          set->size );
        cutilSafeCall(cudaDeviceSynchronize());
        timer->endTimer("fpmf_grad 3rd order");
        break;
      case 4:
        _op_cuda_fpmf_grad_3d<4,num_cells,1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
          (DG_FP *) argGeof.data_d,
          (DG_FP *) arg14.data_d,
          (DG_FP *) arg15.data_d,
          (DG_FP *) arg16.data_d,
          set->size );
        break;
      case 5:
        _op_cuda_fpmf_grad_3d<5,num_cells,1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
          (DG_FP *) argGeof.data_d,
          (DG_FP *) arg14.data_d,
          (DG_FP *) arg15.data_d,
          (DG_FP *) arg16.data_d,
          set->size );
        break;
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
