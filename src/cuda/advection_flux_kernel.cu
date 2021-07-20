//
// auto-generated by op2.py
//

//user function
__device__ void advection_flux_gpu( const double *u, const double *v, double *f0,
                           double *f1, double *f2, double *f3) {
  for(int i = 0; i < 15; i++) {
    f0[i] = u[i] * u[i];
    f1[i] = u[i] * v[i];
    f2[i] = u[i] * v[i];
    f3[i] = v[i] * v[i];
  }

}

// CUDA kernel function
__global__ void op_cuda_advection_flux(
  const double *__restrict arg0,
  const double *__restrict arg1,
  double *arg2,
  double *arg3,
  double *arg4,
  double *arg5,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    advection_flux_gpu(arg0+n*15,
                   arg1+n*15,
                   arg2+n*15,
                   arg3+n*15,
                   arg4+n*15,
                   arg5+n*15);
  }
}


//host stub function
void op_par_loop_advection_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

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
  op_timing_realloc(28);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[28].name      = name;
  OP_kernels[28].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  advection_flux");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_28
      int nthread = OP_BLOCK_SIZE_28;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_advection_flux<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
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
  OP_kernels[28].time     += wall_t2 - wall_t1;
  OP_kernels[28].transfer += (float)set->size * arg0.size;
  OP_kernels[28].transfer += (float)set->size * arg1.size;
  OP_kernels[28].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[28].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[28].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[28].transfer += (float)set->size * arg5.size * 2.0f;
}
