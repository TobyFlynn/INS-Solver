//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void advection_numerical_flux_openacc( const double *fscale, const double *nx,
                                     const double *ny, const double *q0,
                                     const double *q1, double *exQ0,
                                     double *exQ1, double *flux0, double *flux1) {

  double fM[4][3 * 2];
  for(int i = 0; i < 3 * 2; i++) {
    fM[0][i] = q0[FMASK[i]] * q0[FMASK[i]];
    fM[1][i] = q0[FMASK[i]] * q1[FMASK[i]];
    fM[2][i] = q0[FMASK[i]] * q1[FMASK[i]];
    fM[3][i] = q1[FMASK[i]] * q1[FMASK[i]];
  }
  double fP[4][3 * 2];
  for(int i = 0; i < 3 * 2; i++) {
    fP[0][i] = exQ0[i] * exQ0[i];
    fP[1][i] = exQ0[i] * exQ1[i];
    fP[2][i] = exQ0[i] * exQ1[i];
    fP[3][i] = exQ1[i] * exQ1[i];
  }

  double maxVel[3 * 2];
  double max = 0.0;
  for(int i = 0; i < 2; i++) {
    double mVel = q0[FMASK[i]] * nx[i] + q1[FMASK[i]] * ny[i];
    double pVel = exQ0[i] * nx[i] + exQ1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 0; i < 2; i++) {
    maxVel[i] = max;
  }
  max = 0.0;
  for(int i = 2; i < 2 * 2; i++) {
    double mVel = q0[FMASK[i]] * nx[i] + q1[FMASK[i]] * ny[i];
    double pVel = exQ0[i] * nx[i] + exQ1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 2; i < 2 * 2; i++) {
    maxVel[i] = max;
  }
  max = 0.0;
  for(int i = 2 * 2; i < 3 * 2; i++) {
    double mVel = q0[FMASK[i]] * nx[i] + q1[FMASK[i]] * ny[i];
    double pVel = exQ0[i] * nx[i] + exQ1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 2 * 2; i < 3 * 2; i++) {
    maxVel[i] = max;
  }

  for(int i = 0; i < 3 * 2; i++) {
    flux0[i] = 0.5 * fscale[i] * (-nx[i] * (fM[0][i] - fP[0][i]) - ny[i] * (fM[1][i] - fP[1][i]) - maxVel[i] * (exQ0[i] - q0[FMASK[i]]));
    flux1[i] = 0.5 * fscale[i] * (-nx[i] * (fM[2][i] - fP[2][i]) - ny[i] * (fM[3][i] - fP[3][i]) - maxVel[i] * (exQ1[i] - q1[FMASK[i]]));
  }

  for(int i = 0; i < 3 * 2; i++) {
    exQ0[i] = 0.0;
    exQ1[i] = 0.0;
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
  op_timing_realloc(31);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[31].name      = name;
  OP_kernels[31].count    += 1;


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
    double* data7 = (double*)arg7.data_d;
    double* data8 = (double*)arg8.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4,data5,data6,data7,data8)
    for ( int n=0; n<set->size; n++ ){
      advection_numerical_flux_openacc(
        &data0[6*n],
        &data1[6*n],
        &data2[6*n],
        &data3[3*n],
        &data4[3*n],
        &data5[6*n],
        &data6[6*n],
        &data7[6*n],
        &data8[6*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[31].time     += wall_t2 - wall_t1;
  OP_kernels[31].transfer += (float)set->size * arg0.size;
  OP_kernels[31].transfer += (float)set->size * arg1.size;
  OP_kernels[31].transfer += (float)set->size * arg2.size;
  OP_kernels[31].transfer += (float)set->size * arg3.size;
  OP_kernels[31].transfer += (float)set->size * arg4.size;
  OP_kernels[31].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[31].transfer += (float)set->size * arg6.size * 2.0f;
  OP_kernels[31].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[31].transfer += (float)set->size * arg8.size * 2.0f;
}
