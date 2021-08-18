//
// auto-generated by op2.py
//

//user function
__device__ void pressure_mu_gpu( const double *mu, double *curl) {
  for(int i = 0; i < 10; i++) {
    curl[i] *= mu[i];
  }

}

// CUDA kernel function
__global__ void op_cuda_pressure_mu(
  const double *__restrict arg0,
  double *arg1,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    pressure_mu_gpu(arg0+n*10,
                arg1+n*10);
  }
}


//host stub function
void op_par_loop_pressure_mu(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(34);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[34].name      = name;
  OP_kernels[34].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pressure_mu");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_34
      int nthread = OP_BLOCK_SIZE_34;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_pressure_mu<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[34].time     += wall_t2 - wall_t1;
  OP_kernels[34].transfer += (float)set->size * arg0.size;
  OP_kernels[34].transfer += (float)set->size * arg1.size * 2.0f;
}
