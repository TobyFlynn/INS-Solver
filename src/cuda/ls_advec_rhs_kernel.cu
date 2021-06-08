//
// auto-generated by op2.py
//

//user function
__device__ void ls_advec_rhs_gpu( const double *dFdr, const double *dFds,
                         const double *dGdr, const double *dGds,
                         const double *rx, const double *ry, const double *sx,
                         const double *sy, const double *q, double *exQ,
                         const double *u, const double *v, const double *fscale,
                         const double *nx, const double *ny, double *nFlux,
                         double *output) {
  for(int i = 0; i < 15; i++) {
    output[i] = rx[i] * dFdr[i] + sx[i] * dFds[i] + ry[i] * dGdr[i] + sy[i] * dGds[i];
  }

  double mQ[15];
  double mF[15];
  double mG[15];
  for(int i = 0; i < 15; i++) {
    int ind = FMASK_cuda[i];
    mQ[i] = q[ind];
    mF[i] = u[ind] * q[ind];
    mG[i] = v[ind] * q[ind];
  }

  double pF[15];
  double pG[15];
  for(int i = 0; i < 15; i++) {
    int ind = FMASK_cuda[i];
    pF[i]  = u[ind] * exQ[i];
    pG[i]  = v[ind] * exQ[i];
  }

  for(int i = 0; i < 15; i++) {
    int ind = FMASK_cuda[i];



    nFlux[i] = (nx[i] * u[ind] + ny[i] * v[ind]) * (q[ind] + exQ[i]);
    nFlux[i] += fabs(nx[i] * u[ind] + ny[i] * v[ind]) * (q[ind] - exQ[i]);
    nFlux[i] *= 0.5 * fscale[i];
  }

  for(int i = 0; i < 15; i++) {
    exQ[i] = 0.0;
  }

}

// CUDA kernel function
__global__ void op_cuda_ls_advec_rhs(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  const double *__restrict arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  const double *__restrict arg8,
  double *arg9,
  const double *__restrict arg10,
  const double *__restrict arg11,
  const double *__restrict arg12,
  const double *__restrict arg13,
  const double *__restrict arg14,
  double *arg15,
  double *arg16,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    ls_advec_rhs_gpu(arg0+n*15,
                 arg1+n*15,
                 arg2+n*15,
                 arg3+n*15,
                 arg4+n*15,
                 arg5+n*15,
                 arg6+n*15,
                 arg7+n*15,
                 arg8+n*15,
                 arg9+n*15,
                 arg10+n*15,
                 arg11+n*15,
                 arg12+n*15,
                 arg13+n*15,
                 arg14+n*15,
                 arg15+n*15,
                 arg16+n*15);
  }
}


//host stub function
void op_par_loop_ls_advec_rhs(char const *name, op_set set,
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
  op_timing_realloc(55);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[55].name      = name;
  OP_kernels[55].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_advec_rhs");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_55
      int nthread = OP_BLOCK_SIZE_55;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_ls_advec_rhs<<<nblocks,nthread>>>(
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
  OP_kernels[55].time     += wall_t2 - wall_t1;
  OP_kernels[55].transfer += (float)set->size * arg0.size;
  OP_kernels[55].transfer += (float)set->size * arg1.size;
  OP_kernels[55].transfer += (float)set->size * arg2.size;
  OP_kernels[55].transfer += (float)set->size * arg3.size;
  OP_kernels[55].transfer += (float)set->size * arg4.size;
  OP_kernels[55].transfer += (float)set->size * arg5.size;
  OP_kernels[55].transfer += (float)set->size * arg6.size;
  OP_kernels[55].transfer += (float)set->size * arg7.size;
  OP_kernels[55].transfer += (float)set->size * arg8.size;
  OP_kernels[55].transfer += (float)set->size * arg9.size * 2.0f;
  OP_kernels[55].transfer += (float)set->size * arg10.size;
  OP_kernels[55].transfer += (float)set->size * arg11.size;
  OP_kernels[55].transfer += (float)set->size * arg12.size;
  OP_kernels[55].transfer += (float)set->size * arg13.size;
  OP_kernels[55].transfer += (float)set->size * arg14.size;
  OP_kernels[55].transfer += (float)set->size * arg15.size * 2.0f;
  OP_kernels[55].transfer += (float)set->size * arg16.size * 2.0f;
}
