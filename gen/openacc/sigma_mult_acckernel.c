//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void sigma_mult_openacc( const double *eps, double *sigx, double *sigy,
                       double *fx, double *fy, double *diffF) {
  for(int i = 0; i < 6; i++) {
    sigx[i] *= *eps;
    sigy[i] *= *eps;
  }

  for(int i = 0; i < 12; i++) {
    fx[i] = 0.0;
    fy[i] = 0.0;
    diffF[i] = 0.0;
  }
}

// host stub function
void op_par_loop_sigma_mult(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  double*arg0h = (double *)arg0.data;
  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(60);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[60].name      = name;
  OP_kernels[60].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  sigma_mult");
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
    #pragma acc parallel loop independent deviceptr(data1,data2,data3,data4,data5)
    for ( int n=0; n<set->size; n++ ){
      sigma_mult_openacc(
        &arg0_l,
        &data1[6*n],
        &data2[6*n],
        &data3[12*n],
        &data4[12*n],
        &data5[12*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[60].time     += wall_t2 - wall_t1;
  OP_kernels[60].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[60].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[60].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[60].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[60].transfer += (float)set->size * arg5.size * 2.0f;
}
