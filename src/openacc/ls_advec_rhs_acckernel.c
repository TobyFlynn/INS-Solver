//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void ls_advec_rhs_openacc( const double *dFdr, const double *dFds,
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
    int ind = FMASK[i];
    mQ[i] = q[ind];
    mF[i] = u[ind] * q[ind];
    mG[i] = v[ind] * q[ind];
  }

  double pF[15];
  double pG[15];
  for(int i = 0; i < 15; i++) {
    int ind = FMASK[i];
    pF[i]  = u[ind] * exQ[i];
    pG[i]  = v[ind] * exQ[i];
  }

  for(int i = 0; i < 15; i++) {
    int ind = FMASK[i];



    nFlux[i] = (nx[i] * u[ind] + ny[i] * v[ind]) * (q[ind] + exQ[i]);
    nFlux[i] += fabs(nx[i] * u[ind] + ny[i] * v[ind]) * (q[ind] - exQ[i]);
    nFlux[i] *= 0.5 * fscale[i];
  }

  for(int i = 0; i < 15; i++) {
    exQ[i] = 0.0;
  }
}

// host stub function
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
  op_timing_realloc(53);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[53].name      = name;
  OP_kernels[53].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_advec_rhs");
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
    double* data9 = (double*)arg9.data_d;
    double* data10 = (double*)arg10.data_d;
    double* data11 = (double*)arg11.data_d;
    double* data12 = (double*)arg12.data_d;
    double* data13 = (double*)arg13.data_d;
    double* data14 = (double*)arg14.data_d;
    double* data15 = (double*)arg15.data_d;
    double* data16 = (double*)arg16.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16)
    for ( int n=0; n<set->size; n++ ){
      ls_advec_rhs_openacc(
        &data0[15*n],
        &data1[15*n],
        &data2[15*n],
        &data3[15*n],
        &data4[15*n],
        &data5[15*n],
        &data6[15*n],
        &data7[15*n],
        &data8[15*n],
        &data9[15*n],
        &data10[15*n],
        &data11[15*n],
        &data12[15*n],
        &data13[15*n],
        &data14[15*n],
        &data15[15*n],
        &data16[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[53].time     += wall_t2 - wall_t1;
  OP_kernels[53].transfer += (float)set->size * arg0.size;
  OP_kernels[53].transfer += (float)set->size * arg1.size;
  OP_kernels[53].transfer += (float)set->size * arg2.size;
  OP_kernels[53].transfer += (float)set->size * arg3.size;
  OP_kernels[53].transfer += (float)set->size * arg4.size;
  OP_kernels[53].transfer += (float)set->size * arg5.size;
  OP_kernels[53].transfer += (float)set->size * arg6.size;
  OP_kernels[53].transfer += (float)set->size * arg7.size;
  OP_kernels[53].transfer += (float)set->size * arg8.size;
  OP_kernels[53].transfer += (float)set->size * arg9.size * 2.0f;
  OP_kernels[53].transfer += (float)set->size * arg10.size;
  OP_kernels[53].transfer += (float)set->size * arg11.size;
  OP_kernels[53].transfer += (float)set->size * arg12.size;
  OP_kernels[53].transfer += (float)set->size * arg13.size;
  OP_kernels[53].transfer += (float)set->size * arg14.size;
  OP_kernels[53].transfer += (float)set->size * arg15.size * 2.0f;
  OP_kernels[53].transfer += (float)set->size * arg16.size * 2.0f;
}
