//
// auto-generated by op2.py
//

//user function
__device__ void init_gauss_grad2_gpu( const double *nx, const double *ny, const double *Dx0,
                             const double *Dy0, const double *Dx1, const double *Dy1,
                             const double *Dx2, const double *Dy2, double *d0,
                             double *d1, double *d2) {
  for(int m = 0; m < 7; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      d0[ind] = nx[m] * Dx0[ind] + ny[m] * Dy0[ind];
      d1[ind] = nx[m + 7] * Dx1[ind] + ny[m + 7] * Dy1[ind];
      d2[ind] = nx[m + 14] * Dx2[ind] + ny[m + 14] * Dy2[ind];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_init_gauss_grad2(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  const double *__restrict arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  double *arg8,
  double *arg9,
  double *arg10,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    init_gauss_grad2_gpu(arg0+n*21,
                     arg1+n*21,
                     arg2+n*105,
                     arg3+n*105,
                     arg4+n*105,
                     arg5+n*105,
                     arg6+n*105,
                     arg7+n*105,
                     arg8+n*105,
                     arg9+n*105,
                     arg10+n*105);
  }
}


//host stub function
void op_par_loop_init_gauss_grad2(char const *name, op_set set,
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
  op_arg arg10){

  int nargs = 11;
  op_arg args[11];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(20);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[20].name      = name;
  OP_kernels[20].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_gauss_grad2");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_20
      int nthread = OP_BLOCK_SIZE_20;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_init_gauss_grad2<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      (double *) arg6.data_d,
      (double *) arg7.data_d,
      (double *) arg8.data_d,
      (double *) arg9.data_d,
      (double *) arg10.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[20].time     += wall_t2 - wall_t1;
  OP_kernels[20].transfer += (float)set->size * arg0.size;
  OP_kernels[20].transfer += (float)set->size * arg1.size;
  OP_kernels[20].transfer += (float)set->size * arg2.size;
  OP_kernels[20].transfer += (float)set->size * arg3.size;
  OP_kernels[20].transfer += (float)set->size * arg4.size;
  OP_kernels[20].transfer += (float)set->size * arg5.size;
  OP_kernels[20].transfer += (float)set->size * arg6.size;
  OP_kernels[20].transfer += (float)set->size * arg7.size;
  OP_kernels[20].transfer += (float)set->size * arg8.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg9.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg10.size * 2.0f;
}