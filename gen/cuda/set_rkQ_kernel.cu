//
// auto-generated by op2.py
//

//user function
__device__ void set_rkQ_gpu( const int *stage, const double *dt, const double *s,
                    const double *rk0, const double *rk1, double *rkQ) {
  if(*stage == -1) {
    for(int i = 0; i < 10; i++) {
      rkQ[i] = s[i];
    }
  } else if(*stage == 0) {
    for(int i = 0; i < 10; i++) {
      rkQ[i] = s[i] + (*dt) * rk0[i];
    }
  } else {
    for(int i = 0; i < 10; i++) {
      rkQ[i] = s[i] + 0.5 * (*dt) * (rk0[i] / 4.0 + rk1[i] / 4.0);
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_set_rkQ(
  const int *arg0,
  const double *arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  double *arg5,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    set_rkQ_gpu(arg0,
            arg1,
            arg2+n*10,
            arg3+n*10,
            arg4+n*10,
            arg5+n*10);
  }
}


//host stub function
void op_par_loop_set_rkQ(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int*arg0h = (int *)arg0.data;
  double*arg1h = (double *)arg1.data;
  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(49);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[49].name      = name;
  OP_kernels[49].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  set_rkQ");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(int));
    consts_bytes += ROUND_UP(1*sizeof(double));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg0.data   = OP_consts_h + consts_bytes;
    arg0.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg0.data)[d] = arg0h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    arg1.data   = OP_consts_h + consts_bytes;
    arg1.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg1.data)[d] = arg1h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_49
      int nthread = OP_BLOCK_SIZE_49;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_set_rkQ<<<nblocks,nthread>>>(
      (int *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[49].time     += wall_t2 - wall_t1;
  OP_kernels[49].transfer += (float)set->size * arg2.size;
  OP_kernels[49].transfer += (float)set->size * arg3.size;
  OP_kernels[49].transfer += (float)set->size * arg4.size;
  OP_kernels[49].transfer += (float)set->size * arg5.size * 2.0f;
}
