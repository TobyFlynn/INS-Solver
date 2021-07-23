//
// auto-generated by op2.py
//

//user function
__device__ void advection_numerical_flux_gpu( const double *fscale, const double *nx,
                                     const double *ny, const double *q0,
                                     const double *q1, double *exQ0,
                                     double *exQ1, double *flux0, double *flux1) {

  double fM[4][15];
  for(int i = 0; i < 15; i++) {
    fM[0][i] = q0[FMASK_cuda[i]] * q0[FMASK_cuda[i]];
    fM[1][i] = q0[FMASK_cuda[i]] * q1[FMASK_cuda[i]];
    fM[2][i] = q0[FMASK_cuda[i]] * q1[FMASK_cuda[i]];
    fM[3][i] = q1[FMASK_cuda[i]] * q1[FMASK_cuda[i]];
  }
  double fP[4][15];
  for(int i = 0; i < 15; i++) {
    fP[0][i] = exQ0[i] * exQ0[i];
    fP[1][i] = exQ0[i] * exQ1[i];
    fP[2][i] = exQ0[i] * exQ1[i];
    fP[3][i] = exQ1[i] * exQ1[i];
  }

  double maxVel[15];
  double max = 0.0;
  for(int i = 0; i < 5; i++) {
    double mVel = q0[FMASK_cuda[i]] * nx[i] + q1[FMASK_cuda[i]] * ny[i];
    double pVel = exQ0[i] * nx[i] + exQ1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 0; i < 5; i++) {
    maxVel[i] = max;
  }
  max = 0.0;
  for(int i = 5; i < 10; i++) {
    double mVel = q0[FMASK_cuda[i]] * nx[i] + q1[FMASK_cuda[i]] * ny[i];
    double pVel = exQ0[i] * nx[i] + exQ1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 5; i < 10; i++) {
    maxVel[i] = max;
  }
  max = 0.0;
  for(int i = 10; i < 15; i++) {
    double mVel = q0[FMASK_cuda[i]] * nx[i] + q1[FMASK_cuda[i]] * ny[i];
    double pVel = exQ0[i] * nx[i] + exQ1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 10; i < 15; i++) {
    maxVel[i] = max;
  }

  for(int i = 0; i < 15; i++) {
    flux0[i] = 0.5 * fscale[i] * (-nx[i] * (fM[0][i] - fP[0][i]) - ny[i] * (fM[1][i] - fP[1][i]) - maxVel[i] * (exQ0[i] - q0[FMASK_cuda[i]]));
    flux1[i] = 0.5 * fscale[i] * (-nx[i] * (fM[2][i] - fP[2][i]) - ny[i] * (fM[3][i] - fP[3][i]) - maxVel[i] * (exQ1[i] - q1[FMASK_cuda[i]]));
  }

  for(int i = 0; i < 15; i++) {
    exQ0[i] = 0.0;
    exQ1[i] = 0.0;
  }

}

// CUDA kernel function
__global__ void op_cuda_advection_numerical_flux(
  const double *__restrict arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  double *arg5,
  double *arg6,
  double *arg7,
  double *arg8,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    advection_numerical_flux_gpu(arg0+n*15,
                             arg1+n*15,
                             arg2+n*15,
                             arg3+n*15,
                             arg4+n*15,
                             arg5+n*15,
                             arg6+n*15,
                             arg7+n*15,
                             arg8+n*15);
  }
}


//host stub function
void op_par_loop_advection_numerical_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  int nargs = 9;
  op_arg args[9];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(32);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[32].name      = name;
  OP_kernels[32].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  advection_numerical_flux");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_32
      int nthread = OP_BLOCK_SIZE_32;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_advection_numerical_flux<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      (double *) arg6.data_d,
      (double *) arg7.data_d,
      (double *) arg8.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[32].time     += wall_t2 - wall_t1;
  OP_kernels[32].transfer += (float)set->size * arg0.size;
  OP_kernels[32].transfer += (float)set->size * arg1.size;
  OP_kernels[32].transfer += (float)set->size * arg2.size;
  OP_kernels[32].transfer += (float)set->size * arg3.size;
  OP_kernels[32].transfer += (float)set->size * arg4.size;
  OP_kernels[32].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[32].transfer += (float)set->size * arg6.size * 2.0f;
  OP_kernels[32].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[32].transfer += (float)set->size * arg8.size * 2.0f;
}
