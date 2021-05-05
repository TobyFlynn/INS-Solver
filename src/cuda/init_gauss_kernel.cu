//
// auto-generated by op2.py
//

//user function
__device__ void init_gauss_gpu( double *rx, double *sx, double *ry, double *sy,
                       double *nx, double *ny, double *sJ) {

  double J[21];
  for(int i = 0; i < 21; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  for(int i = 0; i < 21; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }


  for(int i = 0; i < 7; i++) {
    nx[i] = -sx[i];
    ny[i] = -sy[i];
  }

  for(int i = 7; i < 14; i++) {
    nx[i] = rx[i] + sx[i];
    ny[i] = ry[i] + sy[i];
  }

  for(int i = 14; i < 21; i++) {
    nx[i] = -rx[i];
    ny[i] = -ry[i];
  }

  for(int i = 0; i < 21; i++) {
    sJ[i] = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / sJ[i];
    ny[i] = ny[i] / sJ[i];
    sJ[i] = sJ[i] * J[i];
  }

}

// CUDA kernel function
__global__ void op_cuda_init_gauss(
  double *arg0,
  double *arg1,
  double *arg2,
  double *arg3,
  double *arg4,
  double *arg5,
  double *arg6,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    init_gauss_gpu(arg0+n*21,
               arg1+n*21,
               arg2+n*21,
               arg3+n*21,
               arg4+n*21,
               arg5+n*21,
               arg6+n*21);
  }
}


//host stub function
void op_par_loop_init_gauss(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6){

  int nargs = 7;
  op_arg args[7];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(5);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[5].name      = name;
  OP_kernels[5].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_gauss");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_5
      int nthread = OP_BLOCK_SIZE_5;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_init_gauss<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      (double *) arg6.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[5].time     += wall_t2 - wall_t1;
  OP_kernels[5].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[5].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[5].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[5].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[5].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[5].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[5].transfer += (float)set->size * arg6.size * 2.0f;
}
