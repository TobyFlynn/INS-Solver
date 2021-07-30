//
// auto-generated by op2.py
//

//user function
__device__ void gauss_op_gpu( const double *tau, const double *sJ,
                     const double *mD0, double *f0_0, double *f0_1, double *f0_2,
                     const double *mD1, double *f1_0, double *f1_1, double *f1_2,
                     const double *mD2, double *f2_0, double *f2_1, double *f2_2,
                     double *pDy0, double *pDy1, double *pDy2) {

  for(int ind = 0; ind < 3 * 3; ind++) {
    int indT = ((ind * 3) % (3 * 3)) + (ind / 3);
    f0_0[ind] = gFInterp0_g_cuda[indT];
    f0_1[ind] = gFInterp0_g_cuda[indT];
    f0_2[ind] = mD0[indT];
  }

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int ind  = m * 3 + n;
      f0_0[ind] = gaussW_g_cuda[n] * sJ[n] * tau[0] * f0_0[ind];
      f0_1[ind] = gaussW_g_cuda[n] * sJ[n] * f0_1[ind];
      f0_2[ind] = gaussW_g_cuda[n] * sJ[n] * f0_2[ind];
    }
  }

  for(int ind = 0; ind < 3 * 3; ind++) {
    int indT = ((ind * 3) % (3 * 7)) + (ind / 3);
    f1_0[ind] = gFInterp1_g_cuda[indT];
    f1_1[ind] = gFInterp1_g_cuda[indT];
    f1_2[ind] = mD1[indT];
  }

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int ind = m * 3 + n;
      f1_0[ind] = gaussW_g_cuda[n] * sJ[n + 3] * tau[1] * f1_0[ind];
      f1_1[ind] = gaussW_g_cuda[n] * sJ[n + 3] * f1_1[ind];
      f1_2[ind] = gaussW_g_cuda[n] * sJ[n + 3] * f1_2[ind];
    }
  }

  for(int ind = 0; ind < 3 * 3; ind++) {
    int indT = ((ind * 3) % (3 * 3)) + (ind / 3);
    f2_0[ind] = gFInterp2_g_cuda[indT];
    f2_1[ind] = gFInterp2_g_cuda[indT];
    f2_2[ind] = mD2[indT];
  }

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int ind = m * 3 + n;
      f2_0[ind] = gaussW_g_cuda[n] * sJ[n + 2 * 3] * tau[2] * f2_0[ind];
      f2_1[ind] = gaussW_g_cuda[n] * sJ[n + 2 * 3] * f2_1[ind];
      f2_2[ind] = gaussW_g_cuda[n] * sJ[n + 2 * 3] * f2_2[ind];
    }
  }

  for(int i = 0; i < 3 * 3; i++) {
    pDy0[i] = 0.0;
    pDy1[i] = 0.0;
    pDy2[i] = 0.0;
  }

}

// CUDA kernel function
__global__ void op_cuda_gauss_op(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  double *arg3,
  double *arg4,
  double *arg5,
  const double *__restrict arg6,
  double *arg7,
  double *arg8,
  double *arg9,
  const double *__restrict arg10,
  double *arg11,
  double *arg12,
  double *arg13,
  double *arg14,
  double *arg15,
  double *arg16,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    gauss_op_gpu(arg0+n*3,
             arg1+n*9,
             arg2+n*9,
             arg3+n*9,
             arg4+n*9,
             arg5+n*9,
             arg6+n*9,
             arg7+n*9,
             arg8+n*9,
             arg9+n*9,
             arg10+n*9,
             arg11+n*9,
             arg12+n*9,
             arg13+n*9,
             arg14+n*9,
             arg15+n*9,
             arg16+n*9);
  }
}


//host stub function
void op_par_loop_gauss_op(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16){

  int nargs = 17;
  op_arg args[17];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  args[10] = arg10;
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;
  args[16] = arg16;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(10);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[10].name      = name;
  OP_kernels[10].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  gauss_op");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_10
      int nthread = OP_BLOCK_SIZE_10;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_gauss_op<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      (double *) arg6.data_d,
      (double *) arg7.data_d,
      (double *) arg8.data_d,
      (double *) arg9.data_d,
      (double *) arg10.data_d,
      (double *) arg11.data_d,
      (double *) arg12.data_d,
      (double *) arg13.data_d,
      (double *) arg14.data_d,
      (double *) arg15.data_d,
      (double *) arg16.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[10].time     += wall_t2 - wall_t1;
  OP_kernels[10].transfer += (float)set->size * arg0.size;
  OP_kernels[10].transfer += (float)set->size * arg1.size;
  OP_kernels[10].transfer += (float)set->size * arg2.size;
  OP_kernels[10].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg6.size;
  OP_kernels[10].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg8.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg9.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg10.size;
  OP_kernels[10].transfer += (float)set->size * arg11.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg12.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg13.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg14.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg15.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg16.size * 2.0f;
}
