//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void gauss_op_openacc( const double *tau, const double *sJ,
                     const double *mD0, double *f0_0, double *f0_1, double *f0_2,
                     const double *mD1, double *f1_0, double *f1_1, double *f1_2,
                     const double *mD2, double *f2_0, double *f2_1, double *f2_2,
                     double *pDy0, double *pDy1, double *pDy2) {

  for(int ind = 0; ind < 3 * 3; ind++) {
    int indT = ((ind * 3) % (3 * 3)) + (ind / 3);
    f0_0[ind] = gFInterp0_g[indT];
    f0_1[ind] = gFInterp0_g[indT];
    f0_2[ind] = mD0[indT];
  }

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int ind  = m * 3 + n;
      f0_0[ind] = gaussW_g[n] * sJ[n] * tau[0] * f0_0[ind];
      f0_1[ind] = gaussW_g[n] * sJ[n] * f0_1[ind];
      f0_2[ind] = gaussW_g[n] * sJ[n] * f0_2[ind];
    }
  }

  for(int ind = 0; ind < 3 * 3; ind++) {
    int indT = ((ind * 3) % (3 * 7)) + (ind / 3);
    f1_0[ind] = gFInterp1_g[indT];
    f1_1[ind] = gFInterp1_g[indT];
    f1_2[ind] = mD1[indT];
  }

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int ind = m * 3 + n;
      f1_0[ind] = gaussW_g[n] * sJ[n + 3] * tau[1] * f1_0[ind];
      f1_1[ind] = gaussW_g[n] * sJ[n + 3] * f1_1[ind];
      f1_2[ind] = gaussW_g[n] * sJ[n + 3] * f1_2[ind];
    }
  }

  for(int ind = 0; ind < 3 * 3; ind++) {
    int indT = ((ind * 3) % (3 * 3)) + (ind / 3);
    f2_0[ind] = gFInterp2_g[indT];
    f2_1[ind] = gFInterp2_g[indT];
    f2_2[ind] = mD2[indT];
  }

  for(int m = 0; m < 3; m++) {
    for(int n = 0; n < 3; n++) {
      int ind = m * 3 + n;
      f2_0[ind] = gaussW_g[n] * sJ[n + 2 * 3] * tau[2] * f2_0[ind];
      f2_1[ind] = gaussW_g[n] * sJ[n + 2 * 3] * f2_1[ind];
      f2_2[ind] = gaussW_g[n] * sJ[n + 2 * 3] * f2_2[ind];
    }
  }

  for(int i = 0; i < 3 * 3; i++) {
    pDy0[i] = 0.0;
    pDy1[i] = 0.0;
    pDy2[i] = 0.0;
  }
}

// host stub function
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
      gauss_op_openacc(
        &data0[3*n],
        &data1[9*n],
        &data2[9*n],
        &data3[9*n],
        &data4[9*n],
        &data5[9*n],
        &data6[9*n],
        &data7[9*n],
        &data8[9*n],
        &data9[9*n],
        &data10[9*n],
        &data11[9*n],
        &data12[9*n],
        &data13[9*n],
        &data14[9*n],
        &data15[9*n],
        &data16[9*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
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
