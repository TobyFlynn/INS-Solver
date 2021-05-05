//
// auto-generated by op2.py
//

//user function
__device__ void advection_intermediate_vel_gpu( const double *a0, const double *a1,
                                       const double *b0, const double *b1,
                                       const double *g0, const double *dt,
                                       const double *q0, const double *q1,
                                       const double *q0Old, const double *q1Old,
                                       const double *N0, const double *N1,
                                       const double *N0Old, const double *N1Old,
                                       double *q0T, double *q1T) {
  for(int i = 0; i < 15; i++) {
    q0T[i] = ((*a0 * q0[i] + *a1 * q0Old[i]) - *dt * (*b0 * N0[i] + *b1 * N0Old[i])) / *g0;
    q1T[i] = ((*a0 * q1[i] + *a1 * q1Old[i]) - *dt * (*b0 * N1[i] + *b1 * N1Old[i])) / *g0;
  }

}

// CUDA kernel function
__global__ void op_cuda_advection_intermediate_vel(
  const double *arg0,
  const double *arg1,
  const double *arg2,
  const double *arg3,
  const double *arg4,
  const double *arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  const double *__restrict arg8,
  const double *__restrict arg9,
  const double *__restrict arg10,
  const double *__restrict arg11,
  const double *__restrict arg12,
  const double *__restrict arg13,
  double *arg14,
  double *arg15,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    advection_intermediate_vel_gpu(arg0,
                               arg1,
                               arg2,
                               arg3,
                               arg4,
                               arg5,
                               arg6+n*15,
                               arg7+n*15,
                               arg8+n*15,
                               arg9+n*15,
                               arg10+n*15,
                               arg11+n*15,
                               arg12+n*15,
                               arg13+n*15,
                               arg14+n*15,
                               arg15+n*15);
  }
}


//host stub function
void op_par_loop_advection_intermediate_vel(char const *name, op_set set,
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
  op_arg arg15){

  double*arg0h = (double *)arg0.data;
  double*arg1h = (double *)arg1.data;
  double*arg2h = (double *)arg2.data;
  double*arg3h = (double *)arg3.data;
  double*arg4h = (double *)arg4.data;
  double*arg5h = (double *)arg5.data;
  int nargs = 16;
  op_arg args[16];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(50);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[50].name      = name;
  OP_kernels[50].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  advection_intermediate_vel");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(double));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg0.data   = OP_consts_h + consts_bytes;
    arg0.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg0.data)[d] = arg0h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    arg1.data   = OP_consts_h + consts_bytes;
    arg1.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg1.data)[d] = arg1h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    arg2.data   = OP_consts_h + consts_bytes;
    arg2.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg2.data)[d] = arg2h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    arg3.data   = OP_consts_h + consts_bytes;
    arg3.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg3.data)[d] = arg3h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    arg4.data   = OP_consts_h + consts_bytes;
    arg4.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg4.data)[d] = arg4h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    arg5.data   = OP_consts_h + consts_bytes;
    arg5.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg5.data)[d] = arg5h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_50
      int nthread = OP_BLOCK_SIZE_50;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_advection_intermediate_vel<<<nblocks,nthread>>>(
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
      (double *) arg11.data_d,
      (double *) arg12.data_d,
      (double *) arg13.data_d,
      (double *) arg14.data_d,
      (double *) arg15.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[50].time     += wall_t2 - wall_t1;
  OP_kernels[50].transfer += (float)set->size * arg6.size;
  OP_kernels[50].transfer += (float)set->size * arg7.size;
  OP_kernels[50].transfer += (float)set->size * arg8.size;
  OP_kernels[50].transfer += (float)set->size * arg9.size;
  OP_kernels[50].transfer += (float)set->size * arg10.size;
  OP_kernels[50].transfer += (float)set->size * arg11.size;
  OP_kernels[50].transfer += (float)set->size * arg12.size;
  OP_kernels[50].transfer += (float)set->size * arg13.size;
  OP_kernels[50].transfer += (float)set->size * arg14.size * 2.0f;
  OP_kernels[50].transfer += (float)set->size * arg15.size * 2.0f;
}
