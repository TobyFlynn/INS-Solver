//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void pressure_update_vel_openacc( const double *factor, const double *rho, const double *dpdx,
                                const double *dpdy, const double *qt0,
                                const double *qt1, double *qtt0, double *qtt1,
                                double *dpdn, double *prBC) {
  for(int i = 0; i < 10; i++) {
    qtt0[i] = qt0[i] - *factor * dpdx[i] / rho[i];
    qtt1[i] = qt1[i] - *factor * dpdy[i] / rho[i];


    dpdn[i] = 0.0;
  }

  for(int i = 0; i < 18; i++) {
    prBC[i] = 0.0;
  }
}

// host stub function
void op_par_loop_pressure_update_vel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  double*arg0h = (double *)arg0.data;
  int nargs = 10;
  op_arg args[10];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(39);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[39].name      = name;
  OP_kernels[39].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pressure_update_vel");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg0_l = arg0h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    double* data6 = (double*)arg6.data_d;
    double* data7 = (double*)arg7.data_d;
    double* data8 = (double*)arg8.data_d;
    double* data9 = (double*)arg9.data_d;
    #pragma acc parallel loop independent deviceptr(data1,data2,data3,data4,data5,data6,data7,data8,data9)
    for ( int n=0; n<set->size; n++ ){
      pressure_update_vel_openacc(
        &arg0_l,
        &data1[10*n],
        &data2[10*n],
        &data3[10*n],
        &data4[10*n],
        &data5[10*n],
        &data6[10*n],
        &data7[10*n],
        &data8[12*n],
        &data9[18*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[39].time     += wall_t2 - wall_t1;
  OP_kernels[39].transfer += (float)set->size * arg1.size;
  OP_kernels[39].transfer += (float)set->size * arg2.size;
  OP_kernels[39].transfer += (float)set->size * arg3.size;
  OP_kernels[39].transfer += (float)set->size * arg4.size;
  OP_kernels[39].transfer += (float)set->size * arg5.size;
  OP_kernels[39].transfer += (float)set->size * arg6.size * 2.0f;
  OP_kernels[39].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[39].transfer += (float)set->size * arg8.size * 2.0f;
  OP_kernels[39].transfer += (float)set->size * arg9.size * 2.0f;
}
