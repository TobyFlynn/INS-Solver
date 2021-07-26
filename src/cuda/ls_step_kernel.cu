//
// auto-generated by op2.py
//

//user function
__device__ void ls_step_gpu( const double *alpha, const double *s, double *step,
                    double *nu, double *rho) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 15; i++) {
    step[i] = tanh(PI * s[i] / *alpha);
    nu[i] = 0.5 * nu0_cuda * (1.0 + step[i]) + 0.5 * nu1_cuda * (1.0 - step[i]);
    rho[i] = 0.5 * rho0_cuda * (1.0 + step[i]) + 0.5 * rho1_cuda * (1.0 - step[i]);
  }

}

// CUDA kernel function
__global__ void op_cuda_ls_step(
  const double *arg0,
  const double *__restrict arg1,
  double *arg2,
  double *arg3,
  double *arg4,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    ls_step_gpu(arg0,
            arg1+n*15,
            arg2+n*15,
            arg3+n*15,
            arg4+n*15);
  }
}


//host stub function
void op_par_loop_ls_step(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  double*arg0h = (double *)arg0.data;
  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(65);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[65].name      = name;
  OP_kernels[65].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_step");
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

    op_cuda_ls_step<<<nblocks,nthread>>>(
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
  OP_kernels[65].time     += wall_t2 - wall_t1;
  OP_kernels[65].transfer += (float)set->size * arg1.size;
  OP_kernels[65].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[65].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[65].transfer += (float)set->size * arg4.size * 2.0f;
}
