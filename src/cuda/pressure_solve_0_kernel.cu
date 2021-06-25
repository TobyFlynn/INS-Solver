//
// auto-generated by op2.py
//

//user function
__device__ void pressure_solve_0_gpu( const double *J, const double *Dx,
                             const double *Dy, const double *rho,
                             const double *u, double *rhs) {
  double tmpX[46*15];
  double tmpY[46*15];

  for(int m = 0; m < 46; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      tmpX[ind] = J[m] * cubW_g_cuda[m] * Dx[ind] / rho[m];
      tmpY[ind] = J[m] * cubW_g_cuda[m] * Dy[ind] / rho[m];
    }
  }

  double op[15*15];
  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      op[c_ind] = 0.0;
      for(int k = 0; k < 46; k++) {

        int b_ind = k * 15 + j;

        int a_ind = k * 15 + i;
        op[c_ind] += Dx[a_ind] * tmpX[b_ind] + Dy[a_ind] * tmpY[b_ind];
      }
    }
  }

  for(int i = 0; i < 15; i++) {
    rhs[i] = 0.0;
    for(int j = 0; j < 15; j++) {
      int op_ind = i * 15 + j;
      rhs[i] += op[op_ind] * u[j];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_pressure_solve_0(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  double *arg5,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    pressure_solve_0_gpu(arg0+n*46,
                     arg1+n*690,
                     arg2+n*690,
                     arg3+n*46,
                     arg4+n*15,
                     arg5+n*15);
  }
}


//host stub function
void op_par_loop_pressure_solve_0(char const *name, op_set set,
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
  op_timing_realloc(27);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[27].name      = name;
  OP_kernels[27].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pressure_solve_0");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_27
      int nthread = OP_BLOCK_SIZE_27;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_pressure_solve_0<<<nblocks,nthread>>>(
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
  OP_kernels[27].time     += wall_t2 - wall_t1;
  OP_kernels[27].transfer += (float)set->size * arg0.size;
  OP_kernels[27].transfer += (float)set->size * arg1.size;
  OP_kernels[27].transfer += (float)set->size * arg2.size;
  OP_kernels[27].transfer += (float)set->size * arg3.size;
  OP_kernels[27].transfer += (float)set->size * arg4.size;
  OP_kernels[27].transfer += (float)set->size * arg5.size * 2.0f;
}