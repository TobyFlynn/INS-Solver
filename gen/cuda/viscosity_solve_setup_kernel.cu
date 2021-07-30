//
// auto-generated by op2.py
//

//user function
__device__ void viscosity_solve_setup_gpu( const double *mu, const double *rho,
                                  const double *mmConst, double *factor,
                                  double *mmFactor) {
  for(int i = 0; i < 6; i++) {
    factor[i] = mu[i];
    mmFactor[i] = *mmConst * rho[i];
  }

}

// CUDA kernel function
__global__ void op_cuda_viscosity_solve_setup(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *arg2,
  double *arg3,
  double *arg4,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    viscosity_solve_setup_gpu(arg0+n*6,
                          arg1+n*6,
                          arg2,
                          arg3+n*6,
                          arg4+n*6);
  }
}


//host stub function
void op_par_loop_viscosity_solve_setup(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  double*arg2h = (double *)arg2.data;
  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(25);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[25].name      = name;
  OP_kernels[25].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  viscosity_solve_setup");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg2.data   = OP_consts_h + consts_bytes;
    arg2.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg2.data)[d] = arg2h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_25
      int nthread = OP_BLOCK_SIZE_25;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_viscosity_solve_setup<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[25].time     += wall_t2 - wall_t1;
  OP_kernels[25].transfer += (float)set->size * arg0.size;
  OP_kernels[25].transfer += (float)set->size * arg1.size;
  OP_kernels[25].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[25].transfer += (float)set->size * arg4.size * 2.0f;
}
