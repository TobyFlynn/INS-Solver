//
// auto-generated by op2.py
//

//user function
__device__ void init_cubature_grad_gpu( double *rx, double *sx, double *ry,  double *sy,
                               double *Dx, double *Dy) {

  double J[46];
  for(int i = 0; i < 46; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  for(int i = 0; i < 46; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  for(int m = 0; m < 46; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      Dx[ind] = rx[m] * cubVDr_cuda[ind] + sx[m] * cubVDs_cuda[ind];
      Dy[ind] = ry[m] * cubVDr_cuda[ind] + sy[m] * cubVDs_cuda[ind];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_init_cubature_grad(
  double *arg0,
  double *arg1,
  double *arg2,
  double *arg3,
  double *arg4,
  double *arg5,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    init_cubature_grad_gpu(arg0+n*46,
                       arg1+n*46,
                       arg2+n*46,
                       arg3+n*46,
                       arg4+n*690,
                       arg5+n*690);
  }
}


//host stub function
void op_par_loop_init_cubature_grad(char const *name, op_set set,
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
  op_timing_realloc(15);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[15].name      = name;
  OP_kernels[15].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_cubature_grad");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_15
      int nthread = OP_BLOCK_SIZE_15;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_init_cubature_grad<<<nblocks,nthread>>>(
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
  OP_kernels[15].time     += wall_t2 - wall_t1;
  OP_kernels[15].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[15].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[15].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[15].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[15].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[15].transfer += (float)set->size * arg5.size * 2.0f;
}
