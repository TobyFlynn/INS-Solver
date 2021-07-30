//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void init_gauss_grad_openacc( double *rx, double *sx, double *ry,  double *sy,
                            double *Dx0, double *Dy0, double *Dx1, double *Dy1,
                            double *Dx2, double *Dy2) {

  double J[18];
  for(int i = 0; i < 18; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  for(int i = 0; i < 18; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int ind = m * 10 + n;
      Dx0[ind] = rx[m] * gF0Dr_g[ind] + sx[m] * gF0Ds_g[ind];
      Dy0[ind] = ry[m] * gF0Dr_g[ind] + sy[m] * gF0Ds_g[ind];
    }
  }

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int ind = m * 10 + n;
      Dx1[ind] = rx[m + 6] * gF1Dr_g[ind] + sx[m + 6] * gF1Ds_g[ind];
      Dy1[ind] = ry[m + 6] * gF1Dr_g[ind] + sy[m + 6] * gF1Ds_g[ind];
    }
  }

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int ind = m * 10 + n;
      Dx2[ind] = rx[m + 2 * 6] * gF2Dr_g[ind] + sx[m + 2 * 6] * gF2Ds_g[ind];
      Dy2[ind] = ry[m + 2 * 6] * gF2Dr_g[ind] + sy[m + 2 * 6] * gF2Ds_g[ind];
    }
  }
}

// host stub function
void op_par_loop_init_gauss_grad(char const *name, op_set set,
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
  op_timing_realloc(6);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[6].name      = name;
  OP_kernels[6].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_gauss_grad");
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
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4,data5,data6,data7,data8,data9)
    for ( int n=0; n<set->size; n++ ){
      init_gauss_grad_openacc(
        &data0[18*n],
        &data1[18*n],
        &data2[18*n],
        &data3[18*n],
        &data4[60*n],
        &data5[60*n],
        &data6[60*n],
        &data7[60*n],
        &data8[60*n],
        &data9[60*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[6].time     += wall_t2 - wall_t1;
  OP_kernels[6].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg6.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg8.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg9.size * 2.0f;
}
