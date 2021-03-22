//
// auto-generated by op2.py
//

//user function
__device__ void poisson_test_error_gpu( const double *x, const double *y,
                               const double *sol, double *err, double *l2) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 15; i++) {
    double x1 = x[i];
    double y1 = y[i];


    double exact = (1.0 - (x[i] * x[i])) * (2.0 * (y[i] * y[i] * y[i]) - 3.0 * (y[i] * y[i]) + 1.0);

    err[i] = fabs(sol[i] - exact);
    *l2 += err[i] * err[i];
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_test_error(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  double *arg3,
  double *arg4,
  int   set_size ) {

  double arg4_l[1];
  for ( int d=0; d<1; d++ ){
    arg4_l[d]=ZERO_double;
  }

  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    poisson_test_error_gpu(arg0+n*15,
                       arg1+n*15,
                       arg2+n*15,
                       arg3+n*15,
                       arg4_l);
  }

  //global reductions

  for ( int d=0; d<1; d++ ){
    op_reduction<OP_INC>(&arg4[d+blockIdx.x*1],arg4_l[d]);
  }
}


//host stub function
void op_par_loop_poisson_test_error(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  double*arg4h = (double *)arg4.data;
  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(35);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[35].name      = name;
  OP_kernels[35].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_test_error");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_35
      int nthread = OP_BLOCK_SIZE_35;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    //transfer global reduction data to GPU
    int maxblocks = nblocks;
    int reduct_bytes = 0;
    int reduct_size  = 0;
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    reduct_size   = MAX(reduct_size,sizeof(double));
    reallocReductArrays(reduct_bytes);
    reduct_bytes = 0;
    arg4.data   = OP_reduct_h + reduct_bytes;
    arg4.data_d = OP_reduct_d + reduct_bytes;
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        ((double *)arg4.data)[d+b*1] = ZERO_double;
      }
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    mvReductArraysToDevice(reduct_bytes);

    int nshared = reduct_size*nthread;
    op_cuda_poisson_test_error<<<nblocks,nthread,nshared>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      set->size );
    //transfer global reduction data back to CPU
    mvReductArraysToHost(reduct_bytes);
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        arg4h[d] = arg4h[d] + ((double *)arg4.data)[d+b*1];
      }
    }
    arg4.data = (char *)arg4h;
    op_mpi_reduce(&arg4,arg4h);
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[35].time     += wall_t2 - wall_t1;
  OP_kernels[35].transfer += (float)set->size * arg0.size;
  OP_kernels[35].transfer += (float)set->size * arg1.size;
  OP_kernels[35].transfer += (float)set->size * arg2.size;
  OP_kernels[35].transfer += (float)set->size * arg3.size * 2.0f;
}
