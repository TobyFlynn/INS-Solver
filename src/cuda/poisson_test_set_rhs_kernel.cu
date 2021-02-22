//
// auto-generated by op2.py
//

//user function
__device__ void poisson_test_set_rhs_gpu( const double *J, const double *ex, double *rhs) {
  for(int i = 0; i < 15; i++) {
    if(ex[i] < 0.0)
      rhs[FMASK_cuda[i]] = 0.0;
  }

  for(int i = 0; i < 15; i++) {
    rhs[i] *= J[i];
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_test_set_rhs(
  const double *__restrict arg0,
  const double *__restrict arg1,
  double *arg2,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    poisson_test_set_rhs_gpu(arg0+n*15,
                         arg1+n*15,
                         arg2+n*15);
  }
}


//host stub function
void op_par_loop_poisson_test_set_rhs(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(31);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[31].name      = name;
  OP_kernels[31].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_test_set_rhs");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_31
      int nthread = OP_BLOCK_SIZE_31;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_poisson_test_set_rhs<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[31].time     += wall_t2 - wall_t1;
  OP_kernels[31].transfer += (float)set->size * arg0.size;
  OP_kernels[31].transfer += (float)set->size * arg1.size;
  OP_kernels[31].transfer += (float)set->size * arg2.size * 2.0f;
}
