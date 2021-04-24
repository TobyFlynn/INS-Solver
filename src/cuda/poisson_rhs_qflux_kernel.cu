//
// auto-generated by op2.py
//

//user function
__device__ void poisson_rhs_qflux_gpu( const double *nx, const double *ny, const double *sJ,
                              const double *tau, const double *gu, double *uFlux,
                              double *qxFlux, double *qyFlux, double *flux) {
  for(int i = 0; i < 21; i++) {

    flux[i] = nx[i] * qxFlux[i] + ny[i] * qyFlux[i] - tau[i / 7] * (gu[i] - uFlux[i]);
    flux[i] *= gaussW_g_cuda[i % 7] * sJ[i];
  }

  for(int i = 0; i < 21; i++) {
    uFlux[i] = 0.0;
    qxFlux[i] = 0.0;
    qyFlux[i] = 0.0;
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_rhs_qflux(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  double *arg5,
  double *arg6,
  double *arg7,
  double *arg8,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    poisson_rhs_qflux_gpu(arg0+n*21,
                      arg1+n*21,
                      arg2+n*21,
                      arg3+n*3,
                      arg4+n*21,
                      arg5+n*21,
                      arg6+n*21,
                      arg7+n*21,
                      arg8+n*21);
  }
}


//host stub function
void op_par_loop_poisson_rhs_qflux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  int nargs = 9;
  op_arg args[9];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(43);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[43].name      = name;
  OP_kernels[43].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_rhs_qflux");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_43
      int nthread = OP_BLOCK_SIZE_43;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_poisson_rhs_qflux<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      (double *) arg6.data_d,
      (double *) arg7.data_d,
      (double *) arg8.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[43].time     += wall_t2 - wall_t1;
  OP_kernels[43].transfer += (float)set->size * arg0.size;
  OP_kernels[43].transfer += (float)set->size * arg1.size;
  OP_kernels[43].transfer += (float)set->size * arg2.size;
  OP_kernels[43].transfer += (float)set->size * arg3.size;
  OP_kernels[43].transfer += (float)set->size * arg4.size;
  OP_kernels[43].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[43].transfer += (float)set->size * arg6.size * 2.0f;
  OP_kernels[43].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[43].transfer += (float)set->size * arg8.size * 2.0f;
}
