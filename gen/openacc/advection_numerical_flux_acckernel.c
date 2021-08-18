//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void advection_numerical_flux_openacc( const double *fscale, const double *nx,
                                     const double *ny, const double *q0,
                                     const double *q1, double *flux0,
                                     double *flux1) {

  double fM[4][3 * 4];
  for(int i = 0; i < 3 * 4; i++) {
    fM[0][i] = q0[FMASK[i]] * q0[FMASK[i]];
    fM[1][i] = q0[FMASK[i]] * q1[FMASK[i]];
    fM[2][i] = q0[FMASK[i]] * q1[FMASK[i]];
    fM[3][i] = q1[FMASK[i]] * q1[FMASK[i]];
  }
  double fP[4][3 * 4];
  for(int i = 0; i < 3 * 4; i++) {
    fP[0][i] = flux0[i] * flux0[i];
    fP[1][i] = flux0[i] * flux1[i];
    fP[2][i] = flux0[i] * flux1[i];
    fP[3][i] = flux1[i] * flux1[i];
  }

  double maxVel[3 * 4];
  double max = 0.0;
  for(int i = 0; i < 4; i++) {
    double mVel = q0[FMASK[i]] * nx[i] + q1[FMASK[i]] * ny[i];
    double pVel = flux0[i] * nx[i] + flux1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 0; i < 4; i++) {
    maxVel[i] = max;
  }
  max = 0.0;
  for(int i = 4; i < 2 * 4; i++) {
    double mVel = q0[FMASK[i]] * nx[i] + q1[FMASK[i]] * ny[i];
    double pVel = flux0[i] * nx[i] + flux1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 4; i < 2 * 4; i++) {
    maxVel[i] = max;
  }
  max = 0.0;
  for(int i = 2 * 4; i < 3 * 4; i++) {
    double mVel = q0[FMASK[i]] * nx[i] + q1[FMASK[i]] * ny[i];
    double pVel = flux0[i] * nx[i] + flux1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 2 * 4; i < 3 * 4; i++) {
    maxVel[i] = max;
  }

  for(int i = 0; i < 3 * 4; i++) {
    flux0[i] = 0.5 * fscale[i] * (-nx[i] * (fM[0][i] - fP[0][i]) - ny[i] * (fM[1][i] - fP[1][i]) - maxVel[i] * (flux0[i] - q0[FMASK[i]]));
    flux1[i] = 0.5 * fscale[i] * (-nx[i] * (fM[2][i] - fP[2][i]) - ny[i] * (fM[3][i] - fP[3][i]) - maxVel[i] * (flux1[i] - q1[FMASK[i]]));
  }
}

// host stub function
void op_par_loop_advection_numerical_flux(char const *name, op_set set,
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
  op_timing_realloc(35);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[35].name      = name;
  OP_kernels[35].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  advection_numerical_flux");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    double* data6 = (double*)arg6.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4,data5,data6)
    for ( int n=0; n<set->size; n++ ){
      advection_numerical_flux_openacc(
        &data0[12*n],
        &data1[12*n],
        &data2[12*n],
        &data3[10*n],
        &data4[10*n],
        &data5[12*n],
        &data6[12*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[35].time     += wall_t2 - wall_t1;
  OP_kernels[35].transfer += (float)set->size * arg0.size;
  OP_kernels[35].transfer += (float)set->size * arg1.size;
  OP_kernels[35].transfer += (float)set->size * arg2.size;
  OP_kernels[35].transfer += (float)set->size * arg3.size;
  OP_kernels[35].transfer += (float)set->size * arg4.size;
  OP_kernels[35].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[35].transfer += (float)set->size * arg6.size * 2.0f;
}