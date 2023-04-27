#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

#include "dg_compiler_defs.h"
#include "dg_global_constants/dg_mat_constants_3d.h"

template<int dg_np>
__device__ void _fpmf_grad_3d_gpu(const int ind, const DG_FP *dr, const DG_FP *ds,
                             const DG_FP *dt, const DG_FP *u, const DG_FP *fact,
                             const DG_FP *rx, const DG_FP *sx, const DG_FP *tx,
                             const DG_FP *ry, const DG_FP *sy, const DG_FP *ty,
                             const DG_FP *rz, const DG_FP *sz, const DG_FP *tz,
                             DG_FP *ux, DG_FP *uy, DG_FP *uz) {
  if(!(ind < dg_np)) return;

  DG_FP tmp_r = 0.0;
  DG_FP tmp_s = 0.0;
  DG_FP tmp_t = 0.0;
  for(int n = 0; n < dg_np; n++) {
    int mat_ind = DG_MAT_IND(ind, n, dg_np, dg_np);
    tmp_r += dr[mat_ind] * u[n];
    tmp_s += ds[mat_ind] * u[n];
    tmp_t += dt[mat_ind] * u[n];
  }
  ux[ind] = fact[ind] * (*rx * tmp_r + *sx * tmp_s + *tx * tmp_t);
  uy[ind] = fact[ind] * (*ry * tmp_r + *sy * tmp_s + *ty * tmp_t);
  uz[ind] = fact[ind] * (*rz * tmp_r + *sz * tmp_s + *tz * tmp_t);
}

// CUDA kernel function
template<int p, int NUM_CELLS>
__global__ void _op_cuda_fpmf_grad_3d(
  const int *__restrict arg0,
  const DG_FP *__restrict arg4,
  const DG_FP *__restrict argFactor,
  const DG_FP *__restrict arg5,
  const DG_FP *__restrict arg6,
  const DG_FP *__restrict arg7,
  const DG_FP *__restrict arg8,
  const DG_FP *__restrict arg9,
  const DG_FP *__restrict arg10,
  const DG_FP *__restrict arg11,
  const DG_FP *__restrict arg12,
  const DG_FP *__restrict arg13,
  DG_FP *arg14,
  DG_FP *arg15,
  DG_FP *arg16,
  int   set_size ) {
  const int np = (p + 1) * (p + 2) * (p + 3) / 6;
  // Load matrices into shared memory
  __shared__ DG_FP dr_shared[DG_NP * DG_NP];
  __shared__ DG_FP ds_shared[DG_NP * DG_NP];
  __shared__ DG_FP dt_shared[DG_NP * DG_NP];
  __shared__ DG_FP u_shared[NUM_CELLS * DG_NP];

  const int start_ind_mat = (p - 1) * DG_NP * DG_NP;
  for(int i = threadIdx.x; i < DG_NP * DG_NP; i += blockDim.x) {
    dr_shared[i] = dg_Dr_kernel[start_ind_mat + i];
    ds_shared[i] = dg_Ds_kernel[start_ind_mat + i];
    dt_shared[i] = dg_Dt_kernel[start_ind_mat + i];
  }

  __syncthreads();

  //process set elements
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
        u_shared[i] = arg4[start_ind + i];
      }
      // u_shared[threadIdx.x] = arg4[cell_id * DG_NP + node_id];
      __syncthreads();

      _fpmf_grad_3d_gpu<np>(node_id,
                dr_shared,
                ds_shared,
                dt_shared,
                u_shared + local_cell_id * DG_NP, //arg4 + cell_id * DG_NP,
                argFactor + cell_id * DG_NP,
                arg5 + cell_id * 1,
                arg6 + cell_id * 1,
                arg7 + cell_id * 1,
                arg8 + cell_id * 1,
                arg9 + cell_id * 1,
                arg10 + cell_id * 1,
                arg11 + cell_id * 1,
                arg12 + cell_id * 1,
                arg13 + cell_id * 1,
                arg14 + cell_id * DG_NP,
                arg15 + cell_id * DG_NP,
                arg16 + cell_id * DG_NP);
    } else {
      _fpmf_grad_3d_gpu<np>(node_id,
                dr_shared,
                ds_shared,
                dt_shared,
                arg4 + cell_id * DG_NP,
                argFactor + cell_id * DG_NP,
                arg5 + cell_id * 1,
                arg6 + cell_id * 1,
                arg7 + cell_id * 1,
                arg8 + cell_id * 1,
                arg9 + cell_id * 1,
                arg10 + cell_id * 1,
                arg11 + cell_id * 1,
                arg12 + cell_id * 1,
                arg13 + cell_id * 1,
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

  int nargs = 15;
  op_arg args[15];

  args[0] = arg0;
  args[1] = arg4;
  args[2] = argFactor;
  args[3] = arg5;
  args[4] = arg6;
  args[5] = arg7;
  args[6] = arg8;
  args[7] = arg9;
  args[8] = arg10;
  args[9] = arg11;
  args[10] = arg12;
  args[11] = arg13;
  args[12] = arg14;
  args[13] = arg15;
  args[14] = arg16;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  grad_3d");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {
    //set CUDA execution parameters
    const int nthread = 256;
    const int nblocks = 200 < (set->size * DG_NP) / nthread + 1 ? 200 : (set->size * DG_NP) / nthread + 1;
    const int num_cells = (nthread / DG_NP) + 2;

    switch(order) {
      case 1:
        _op_cuda_fpmf_grad_3d<1,num_cells><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
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
        break;
      case 2:
        _op_cuda_fpmf_grad_3d<2,num_cells><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
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
        break;
      case 3:
        timer->startTimer("fpmf_grad 3rd order");
        _op_cuda_fpmf_grad_3d<3,num_cells><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
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
        cutilSafeCall(cudaDeviceSynchronize());
        timer->endTimer("fpmf_grad 3rd order");
        break;
      case 4:
        _op_cuda_fpmf_grad_3d<4,num_cells><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
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
        break;
      case 5:
        _op_cuda_fpmf_grad_3d<5,num_cells><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (DG_FP *) arg4.data_d,
          (DG_FP *) argFactor.data_d,
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
        break;
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
