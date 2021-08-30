//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void pressure_rhs_openacc( const double *b0, const double *b1, const double *g0,
                         const double *dt, const double *J, const double *sJ,
                         const double *dPdN, double *dPdNOld, double *divVelT) {

  double factor = 1.0 / (*dt);
  for(int i = 0; i < 10; i++) {
    divVelT[i] = J[i] * (-divVelT[i] * factor);
  }

  for(int i = 0; i < 3 * 4; i++) {
    dPdNOld[i] = sJ[i] * ((*b0) * dPdN[i] + (*b1) * dPdNOld[i]);
  }
}

// host stub function
void op_par_loop_pressure_rhs(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  double*arg0h = (double *)arg0.data;
  double*arg1h = (double *)arg1.data;
  double*arg2h = (double *)arg2.data;
  double*arg3h = (double *)arg3.data;
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
  op_timing_realloc(39);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[39].name      = name;
  OP_kernels[39].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pressure_rhs");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg0_l = arg0h[0];
  double arg1_l = arg1h[0];
  double arg2_l = arg2h[0];
  double arg3_l = arg3h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    double* data6 = (double*)arg6.data_d;
    double* data7 = (double*)arg7.data_d;
    double* data8 = (double*)arg8.data_d;
    #pragma acc parallel loop independent deviceptr(data4,data5,data6,data7,data8)
    for ( int n=0; n<set->size; n++ ){
      pressure_rhs_openacc(
        &arg0_l,
        &arg1_l,
        &arg2_l,
        &arg3_l,
        &data4[10*n],
        &data5[12*n],
        &data6[12*n],
        &data7[12*n],
        &data8[10*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[39].time     += wall_t2 - wall_t1;
  OP_kernels[39].transfer += (float)set->size * arg4.size;
  OP_kernels[39].transfer += (float)set->size * arg5.size;
  OP_kernels[39].transfer += (float)set->size * arg6.size;
  OP_kernels[39].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[39].transfer += (float)set->size * arg8.size * 2.0f;
}
