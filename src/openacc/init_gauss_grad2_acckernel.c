//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void init_gauss_grad2_openacc( const double *nx, const double *ny, const double *Dx0,
                             const double *Dy0, const double *Dx1, const double *Dy1,
                             const double *Dx2, const double *Dy2, double *d0,
                             double *d1, double *d2) {
  for(int m = 0; m < 7; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      d0[ind] = nx[m] * Dx0[ind] + ny[m] * Dy0[ind];
      d1[ind] = nx[m + 7] * Dx1[ind] + ny[m + 7] * Dy1[ind];
      d2[ind] = nx[m + 14] * Dx2[ind] + ny[m + 14] * Dy2[ind];
    }
  }
}

// host stub function
void op_par_loop_init_gauss_grad2(char const *name, op_set set,
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
  op_arg arg10){

  int nargs = 11;
  op_arg args[11];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(22);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[22].name      = name;
  OP_kernels[22].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_gauss_grad2");
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
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4,data5,data6,data7,data8,data9,data10)
    for ( int n=0; n<set->size; n++ ){
      init_gauss_grad2_openacc(
        &data0[21*n],
        &data1[21*n],
        &data2[105*n],
        &data3[105*n],
        &data4[105*n],
        &data5[105*n],
        &data6[105*n],
        &data7[105*n],
        &data8[105*n],
        &data9[105*n],
        &data10[105*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[22].time     += wall_t2 - wall_t1;
  OP_kernels[22].transfer += (float)set->size * arg0.size;
  OP_kernels[22].transfer += (float)set->size * arg1.size;
  OP_kernels[22].transfer += (float)set->size * arg2.size;
  OP_kernels[22].transfer += (float)set->size * arg3.size;
  OP_kernels[22].transfer += (float)set->size * arg4.size;
  OP_kernels[22].transfer += (float)set->size * arg5.size;
  OP_kernels[22].transfer += (float)set->size * arg6.size;
  OP_kernels[22].transfer += (float)set->size * arg7.size;
  OP_kernels[22].transfer += (float)set->size * arg8.size * 2.0f;
  OP_kernels[22].transfer += (float)set->size * arg9.size * 2.0f;
  OP_kernels[22].transfer += (float)set->size * arg10.size * 2.0f;
}
