//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void init_gauss_openacc( double *rx, double *sx, double *ry, double *sy,
                       double *nx, double *ny, double *sJ) {

  double J[21];
  for(int i = 0; i < 21; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  for(int i = 0; i < 21; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }


  for(int i = 0; i < 7; i++) {
    nx[i] = -sx[i];
    ny[i] = -sy[i];
  }

  for(int i = 7; i < 14; i++) {
    nx[i] = rx[i] + sx[i];
    ny[i] = ry[i] + sy[i];
  }

  for(int i = 14; i < 21; i++) {
    nx[i] = -rx[i];
    ny[i] = -ry[i];
  }

  for(int i = 0; i < 21; i++) {
    sJ[i] = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / sJ[i];
    ny[i] = ny[i] / sJ[i];
    sJ[i] = sJ[i] * J[i];
  }
}

// host stub function
void op_par_loop_init_gauss(char const *name, op_set set,
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
  op_timing_realloc(20);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[20].name      = name;
  OP_kernels[20].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_gauss");
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
      init_gauss_openacc(
        &data0[21*n],
        &data1[21*n],
        &data2[21*n],
        &data3[21*n],
        &data4[21*n],
        &data5[21*n],
        &data6[21*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[20].time     += wall_t2 - wall_t1;
  OP_kernels[20].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[20].transfer += (float)set->size * arg6.size * 2.0f;
}
