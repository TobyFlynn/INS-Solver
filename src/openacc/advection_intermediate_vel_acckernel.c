//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void advection_intermediate_vel_openacc( const double *a0, const double *a1,
                                       const double *b0, const double *b1,
                                       const double *g0, const double *dt,
                                       const double *q0, const double *q1,
                                       const double *q0Old, const double *q1Old,
                                       const double *N0, const double *N1,
                                       const double *N0Old, const double *N1Old,
                                       double *q0T, double *q1T) {
  for(int i = 0; i < 15; i++) {
    q0T[i] = ((*a0 * q0[i] + *a1 * q0Old[i]) - *dt * (*b0 * N0[i] + *b1 * N0Old[i])) / *g0;
    q1T[i] = ((*a0 * q1[i] + *a1 * q1Old[i]) - *dt * (*b0 * N1[i] + *b1 * N1Old[i])) / *g0;
  }
}

// host stub function
void op_par_loop_advection_intermediate_vel(char const *name, op_set set,
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
  op_arg arg15){

  double*arg0h = (double *)arg0.data;
  double*arg1h = (double *)arg1.data;
  double*arg2h = (double *)arg2.data;
  double*arg3h = (double *)arg3.data;
  double*arg4h = (double *)arg4.data;
  double*arg5h = (double *)arg5.data;
  int nargs = 16;
  op_arg args[16];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(45);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[45].name      = name;
  OP_kernels[45].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  advection_intermediate_vel");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg0_l = arg0h[0];
  double arg1_l = arg1h[0];
  double arg2_l = arg2h[0];
  double arg3_l = arg3h[0];
  double arg4_l = arg4h[0];
  double arg5_l = arg5h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

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
    #pragma acc parallel loop independent deviceptr(data6,data7,data8,data9,data10,data11,data12,data13,data14,data15)
    for ( int n=0; n<set->size; n++ ){
      advection_intermediate_vel_openacc(
        &arg0_l,
        &arg1_l,
        &arg2_l,
        &arg3_l,
        &arg4_l,
        &arg5_l,
        &data6[15*n],
        &data7[15*n],
        &data8[15*n],
        &data9[15*n],
        &data10[15*n],
        &data11[15*n],
        &data12[15*n],
        &data13[15*n],
        &data14[15*n],
        &data15[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[45].time     += wall_t2 - wall_t1;
  OP_kernels[45].transfer += (float)set->size * arg6.size;
  OP_kernels[45].transfer += (float)set->size * arg7.size;
  OP_kernels[45].transfer += (float)set->size * arg8.size;
  OP_kernels[45].transfer += (float)set->size * arg9.size;
  OP_kernels[45].transfer += (float)set->size * arg10.size;
  OP_kernels[45].transfer += (float)set->size * arg11.size;
  OP_kernels[45].transfer += (float)set->size * arg12.size;
  OP_kernels[45].transfer += (float)set->size * arg13.size;
  OP_kernels[45].transfer += (float)set->size * arg14.size * 2.0f;
  OP_kernels[45].transfer += (float)set->size * arg15.size * 2.0f;
}
