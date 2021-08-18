//
// auto-generated by op2.py
//

//user function
__device__ void init_gauss_grad_neighbour_gpu( const int *reverse, double *rx,
                            double *sx, double *ry,  double *sy,
                            double *Dx0, double *Dy0, double *Dx1, double *Dy1,
                            double *Dx2, double *Dy2) {

  double J[18];
  for(int i = 0; i < 18; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  for(int i = 0; i < 18; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  if(reverse[0]) {
    for(int m = 0; m < 6; m++) {
      for(int n = 0; n < 10; n++) {
        int ind = m * 10 + n;
        Dx0[ind] = rx[m] * gF0DrR_g_cuda[ind] + sx[m] * gF0DsR_g_cuda[ind];
        Dy0[ind] = ry[m] * gF0DrR_g_cuda[ind] + sy[m] * gF0DsR_g_cuda[ind];
      }
    }
  } else {
    for(int m = 0; m < 6; m++) {
      for(int n = 0; n < 10; n++) {
        int ind = m * 10 + n;
        Dx0[ind] = rx[m] * gF0Dr_g_cuda[ind] + sx[m] * gF0Ds_g_cuda[ind];
        Dy0[ind] = ry[m] * gF0Dr_g_cuda[ind] + sy[m] * gF0Ds_g_cuda[ind];
      }
    }
  }

  if(reverse[1]) {
    for(int m = 0; m < 6; m++) {
      for(int n = 0; n < 10; n++) {
        int ind = m * 10 + n;
        Dx1[ind] = rx[m + 6] * gF1DrR_g_cuda[ind] + sx[m + 6] * gF1DsR_g_cuda[ind];
        Dy1[ind] = ry[m + 6] * gF1DrR_g_cuda[ind] + sy[m + 6] * gF1DsR_g_cuda[ind];
      }
    }
  } else {
    for(int m = 0; m < 6; m++) {
      for(int n = 0; n < 10; n++) {
        int ind = m * 10 + n;
        Dx1[ind] = rx[m + 6] * gF1Dr_g_cuda[ind] + sx[m + 6] * gF1Ds_g_cuda[ind];
        Dy1[ind] = ry[m + 6] * gF1Dr_g_cuda[ind] + sy[m + 6] * gF1Ds_g_cuda[ind];
      }
    }
  }

  if(reverse[2]) {
    for(int m = 0; m < 6; m++) {
      for(int n = 0; n < 10; n++) {
        int ind = m * 10 + n;
        Dx2[ind] = rx[m + 2 * 6] * gF2DrR_g_cuda[ind] + sx[m + 2 * 6] * gF2DsR_g_cuda[ind];
        Dy2[ind] = ry[m + 2 * 6] * gF2DrR_g_cuda[ind] + sy[m + 2 * 6] * gF2DsR_g_cuda[ind];
      }
    }
  } else {
    for(int m = 0; m < 6; m++) {
      for(int n = 0; n < 10; n++) {
        int ind = m * 10 + n;
        Dx2[ind] = rx[m + 2 * 6] * gF2Dr_g_cuda[ind] + sx[m + 2 * 6] * gF2Ds_g_cuda[ind];
        Dy2[ind] = ry[m + 2 * 6] * gF2Dr_g_cuda[ind] + sy[m + 2 * 6] * gF2Ds_g_cuda[ind];
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_init_gauss_grad_neighbour(
  const int *__restrict arg0,
  double *arg1,
  double *arg2,
  double *arg3,
  double *arg4,
  double *arg5,
  double *arg6,
  double *arg7,
  double *arg8,
  double *arg9,
  double *arg10,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    init_gauss_grad_neighbour_gpu(arg0+n*3,
                              arg1+n*18,
                              arg2+n*18,
                              arg3+n*18,
                              arg4+n*18,
                              arg5+n*60,
                              arg6+n*60,
                              arg7+n*60,
                              arg8+n*60,
                              arg9+n*60,
                              arg10+n*60);
  }
}


//host stub function
void op_par_loop_init_gauss_grad_neighbour(char const *name, op_set set,
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
  op_arg arg10){

  int nargs = 11;
  op_arg args[11];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(10);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[10].name      = name;
  OP_kernels[10].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_gauss_grad_neighbour");
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

    op_cuda_init_gauss_grad_neighbour<<<nblocks,nthread>>>(
      (int *) arg0.data_d,
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
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[10].time     += wall_t2 - wall_t1;
  OP_kernels[10].transfer += (float)set->size * arg0.size;
  OP_kernels[10].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg6.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg8.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg9.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg10.size * 2.0f;
}