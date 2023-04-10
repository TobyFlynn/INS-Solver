#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

#include "dg_compiler_defs.h"

__device__ void _p_multigrid_relaxation_chebyshev_0_gpu(const int node, const double *factor,
                                const double *b,
                                const double *diag, double *res, double *out) {
  res[node] = (b[node] - res[node]) / diag[node];
  out[node] = *factor * res[node];
}

template<int p>
__global__ void _op_cuda_p_multigrid_relaxation_chebyshev_0(
  const double *arg0,
  const int *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  double *arg4,
  double *arg5,
  int   set_size ) {

  const int np = (p + 1) * (p + 2) * (p + 3) / 6;

  //process set elements
  for (int n = threadIdx.x + blockIdx.x * blockDim.x;
       n < set_size * np; n += blockDim.x * gridDim.x){

    //user-supplied kernel call
    const int node = n % np;
    const int cell = n / np;
    _p_multigrid_relaxation_chebyshev_0_gpu(node, arg0,
                                       arg2+cell*DG_NP,
                                       arg3+cell*DG_NP,
                                       arg4+cell*DG_NP,
                                       arg5+cell*DG_NP);
  }
}


//host stub function
void custom_kernel_p_multigrid_relaxation_chebyshev_0(const int order, char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  double*arg0h = (double *)arg0.data;
  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  p_multigrid_relaxation_chebyshev_0");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg0.data   = OP_consts_h + consts_bytes;
    arg0.data_d = OP_consts_d + consts_bytes;
    memcpy(arg0.data, arg0h, 1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(double));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    const int nthread = 256;
    const int nblocks = 200 < (set->size * DG_NP) / nthread + 1 ? 200 : (set->size * DG_NP) / nthread + 1;

    switch(order) {
      case 1:
        _op_cuda_p_multigrid_relaxation_chebyshev_0<1><<<nblocks,nthread>>>(
          (double *) arg0.data_d,
          (int *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          set->size );
        break;
      case 2:
        _op_cuda_p_multigrid_relaxation_chebyshev_0<2><<<nblocks,nthread>>>(
          (double *) arg0.data_d,
          (int *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          set->size );
        break;
      case 3:
        _op_cuda_p_multigrid_relaxation_chebyshev_0<3><<<nblocks,nthread>>>(
          (double *) arg0.data_d,
          (int *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          set->size );
        break;
      case 4:
        _op_cuda_p_multigrid_relaxation_chebyshev_0<4><<<nblocks,nthread>>>(
          (double *) arg0.data_d,
          (int *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          set->size );
        break;
      case 5:
        _op_cuda_p_multigrid_relaxation_chebyshev_0<5><<<nblocks,nthread>>>(
          (double *) arg0.data_d,
          (int *) arg1.data_d,
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
