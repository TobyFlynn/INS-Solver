#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

#include "dg_compiler_defs.h"

template<int p, int dg_np>
__device__ void _fpmf_3d_mult_mm_gpu(const int node, const double *mass,
                                    const double *J, const double *mm_factor,
                                    const double *in, double *out) {
  const double *mass_mat = &mass[(p - 1) * 20 * 20];
  DG_FP tmp = 0.0;
  for(int n = 0; n < dg_np; n++) {
    int ind = DG_MAT_IND(node, n, dg_np, dg_np);
    tmp += mm_factor[n] * mass_mat[ind] * in[n];
  }
  out[node] += tmp * J[0];
}

// CUDA kernel function
template<int p>
__global__ void _op_cuda_fpmf_3d_mult_mm(
  const int *__restrict arg0,
  const double *arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  double *arg5,
  int   set_size ) {

  const int np = (p + 1) * (p + 2) * (p + 3) / 6;
  //process set elements
  for(int n = threadIdx.x + blockIdx.x * blockDim.x;
      n < set_size * np; n += blockDim.x * gridDim.x){

    //user-supplied kernel call
    const int node = n % np;
    const int cell = n / np;
    _fpmf_3d_mult_mm_gpu<p,np>(node,
                    arg1,
                    arg2+cell*1,
                    arg3+cell*DG_NP,
                    arg4+cell*DG_NP,
                    arg5+cell*DG_NP);
  }
}


//host stub function
void custom_kernel_fpmf_3d_mult_mm(const int order, char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  double*arg1h = (double *)arg1.data;
  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  fpmf_3d_mult_mm");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1200*sizeof(double));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg1.data   = OP_consts_h + consts_bytes;
    arg1.data_d = OP_consts_d + consts_bytes;
    memcpy(arg1.data, arg1h, 1200*sizeof(double));
    consts_bytes += ROUND_UP(1200*sizeof(double));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    const int nthread = 256;
    const int nblocks = 200 < (set->size * DG_NP) / nthread + 1 ? 200 : (set->size * DG_NP) / nthread + 1;

    switch(order) {
      case 1:
        _op_cuda_fpmf_3d_mult_mm<1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          set->size );
        break;
      case 2:
        _op_cuda_fpmf_3d_mult_mm<2><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          set->size );
        break;
      case 3:
        _op_cuda_fpmf_3d_mult_mm<3><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          set->size );
        break;
      case 4:
        _op_cuda_fpmf_3d_mult_mm<4><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          set->size );
        break;
      case 5:
        _op_cuda_fpmf_3d_mult_mm<5><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          set->size );
        break;
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
