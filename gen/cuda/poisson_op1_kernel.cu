//
// auto-generated by op2.py
//

//user function
__device__ void poisson_op1_gpu( const double *J, const double *Dx, const double *Dy,
                        const double *factor, double *op) {
  double tmpX[36 * 10];
  double tmpY[36 * 10];

  for(int m = 0; m < 36; m++) {
    for(int n = 0; n < 10; n++) {
      int ind = m * 10 + n;
      tmpX[ind] = J[m] * cubW_g_cuda[m] * Dx[ind] * factor[m];
      tmpY[ind] = J[m] * cubW_g_cuda[m] * Dy[ind] * factor[m];
    }
  }

  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      op[c_ind] = 0.0;
      for(int k = 0; k < 36; k++) {

        int b_ind = k * 10 + j;

        int a_ind = k * 10 + i;
        op[c_ind] += Dx[a_ind] * tmpX[b_ind] + Dy[a_ind] * tmpY[b_ind];
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_op1(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  double *arg4,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    poisson_op1_gpu(arg0+n*36,
                arg1+n*360,
                arg2+n*360,
                arg3+n*36,
                arg4+n*100);
  }
}


//host stub function
void op_par_loop_poisson_op1(char const *name, op_set set,
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
  op_timing_realloc(22);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[22].name      = name;
  OP_kernels[22].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_op1");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_22
      int nthread = OP_BLOCK_SIZE_22;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_poisson_op1<<<nblocks,nthread>>>(
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
  OP_kernels[22].time     += wall_t2 - wall_t1;
  OP_kernels[22].transfer += (float)set->size * arg0.size;
  OP_kernels[22].transfer += (float)set->size * arg1.size;
  OP_kernels[22].transfer += (float)set->size * arg2.size;
  OP_kernels[22].transfer += (float)set->size * arg3.size;
  OP_kernels[22].transfer += (float)set->size * arg4.size * 2.0f;
}