//
// auto-generated by op2.py
//

//user function
__device__ void sigma_mult_gpu( const double *eps, double *sigx, double *sigy,
                       double *diffF) {
  for(int i = 0; i < 10; i++) {
    sigx[i] *= *eps;
    sigy[i] *= *eps;
  }

  for(int i = 0; i < 18; i++) {
    diffF[i] = 0.0;
  }

}

// CUDA kernel function
__global__ void op_cuda_sigma_mult(
  const double *arg0,
  double *arg1,
  double *arg2,
  double *arg3,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    sigma_mult_gpu(arg0,
               arg1+n*10,
               arg2+n*10,
               arg3+n*18);
  }
}


//host stub function
void op_par_loop_sigma_mult(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  double*arg0h = (double *)arg0.data;
  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(65);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[65].name      = name;
  OP_kernels[65].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  sigma_mult");
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
    for ( int d=0; d<1; d++ ){
      ((double *)arg0.data)[d] = arg0h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_65
      int nthread = OP_BLOCK_SIZE_65;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_sigma_mult<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[65].time     += wall_t2 - wall_t1;
  OP_kernels[65].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[65].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[65].transfer += (float)set->size * arg3.size * 2.0f;
}
