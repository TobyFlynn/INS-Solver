#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

#include "dg_compiler_defs.h"

__device__ void _p_multigrid_relaxation_chebyshev_1_gpu(const int node, const double *in,
                                               double *out) {
  out[node] += in[node];
}

// CUDA kernel function
template<int p>
__global__ void _op_cuda_p_multigrid_relaxation_chebyshev_1(
  const int *__restrict arg0,
  const double *__restrict arg1,
  double *arg2,
  int   set_size ) {

  const int np = (p + 1) * (p + 2) * (p + 3) / 6;

  //process set elements
  for(int n = threadIdx.x + blockIdx.x * blockDim.x;
      n < set_size * np; n += blockDim.x * gridDim.x){

    const int node = n % np;
    const int cell = n / np;
    _p_multigrid_relaxation_chebyshev_1_gpu(node, arg1+cell*DG_NP, arg2+cell*DG_NP);
  }
}

void custom_kernel_p_multigrid_relaxation_chebyshev_1(const int order, char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  p_multigrid_relaxation_chebyshev_1");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2, 0);
  if (set_size > 0) {

    //set CUDA execution parameters
    const int nthread = 256;
    const int nblocks = 200 < (set->size * DG_NP) / nthread + 1 ? 200 : (set->size * DG_NP) / nthread + 1;

    switch(order) {
      case 1:
        _op_cuda_p_multigrid_relaxation_chebyshev_1<1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          set->size );
        break;
      case 2:
        _op_cuda_p_multigrid_relaxation_chebyshev_1<2><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          set->size );
        break;
      case 3:
        _op_cuda_p_multigrid_relaxation_chebyshev_1<3><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          set->size );
        break;
      case 4:
        _op_cuda_p_multigrid_relaxation_chebyshev_1<4><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          set->size );
        break;
      case 5:
        _op_cuda_p_multigrid_relaxation_chebyshev_1<5><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          set->size );
        break;
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
