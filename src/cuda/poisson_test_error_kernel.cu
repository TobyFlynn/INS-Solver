//
// auto-generated by op2.py
//

//user function
__device__ void poisson_test_error_gpu( const double *x, const double *y,
                               const double *sol, double *err) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 15; i++) {
    double x1 = x[i];
    double y1 = y[i];



    double exact = y1 * (1.0 - y1) * x1 * x1 * x1;
    err[i] = fabs(sol[i] - exact);

  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_test_error(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  double *arg3,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    poisson_test_error_gpu(arg0+n*15,
                       arg1+n*15,
                       arg2+n*15,
                       arg3+n*15);
  }
}


//host stub function
void op_par_loop_poisson_test_error(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(50);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[50].name      = name;
  OP_kernels[50].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_test_error");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_50
      int nthread = OP_BLOCK_SIZE_50;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_poisson_test_error<<<nblocks,nthread>>>(
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
  OP_kernels[50].time     += wall_t2 - wall_t1;
  OP_kernels[50].transfer += (float)set->size * arg0.size;
  OP_kernels[50].transfer += (float)set->size * arg1.size;
  OP_kernels[50].transfer += (float)set->size * arg2.size;
  OP_kernels[50].transfer += (float)set->size * arg3.size * 2.0f;
}
