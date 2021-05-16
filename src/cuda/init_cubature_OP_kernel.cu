//
// auto-generated by op2.py
//

//user function
__device__ void init_cubature_OP_gpu( const double *J, const double *Dx, const double *Dy,
                             double *temp, double *temp2) {
  for(int m = 0; m < 46; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      temp[ind] = J[m] * cubW_g_cuda[m] * Dx[ind];
      temp2[ind] = J[m] * cubW_g_cuda[m] * Dy[ind];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_init_cubature_OP(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  double *arg3,
  double *arg4,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    init_cubature_OP_gpu(arg0+n*46,
                     arg1+n*690,
                     arg2+n*690,
                     arg3+n*690,
                     arg4+n*690);
  }
}


//host stub function
void op_par_loop_init_cubature_OP(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(4);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[4].name      = name;
  OP_kernels[4].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_cubature_OP");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_4
      int nthread = OP_BLOCK_SIZE_4;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_init_cubature_OP<<<nblocks,nthread>>>(
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
  OP_kernels[4].time     += wall_t2 - wall_t1;
  OP_kernels[4].transfer += (float)set->size * arg0.size;
  OP_kernels[4].transfer += (float)set->size * arg1.size;
  OP_kernels[4].transfer += (float)set->size * arg2.size;
  OP_kernels[4].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg4.size * 2.0f;
}
