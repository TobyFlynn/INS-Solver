//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void set_ic_openacc( const int *problem, const double *x, const double *y,
                   const double *nu, double *q0, double *q1) {
  const double PI = 3.141592653589793238463;
  if(*problem == 0) {
    for(int i = 0; i < 15; i++) {
      q0[i] = 0.0;
      q1[i] = 0.0;
    }
  } else {
    for(int i = 0; i < 15; i++) {
      q0[i] = -sin(2.0 * PI * y[i]) * exp(-nu[i] * 4.0 * PI * PI * 0.0);
      q1[i] = sin(2.0 * PI * x[i]) * exp(-nu[i] * 4.0 * PI * PI * 0.0);
    }
  }
}

// host stub function
void op_par_loop_set_ic(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int*arg0h = (int *)arg0.data;
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
  op_timing_realloc(49);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[49].name      = name;
  OP_kernels[49].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  set_ic");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  int arg0_l = arg0h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    #pragma acc parallel loop independent deviceptr(data1,data2,data3,data4,data5)
    for ( int n=0; n<set->size; n++ ){
      set_ic_openacc(
        &arg0_l,
        &data1[15*n],
        &data2[15*n],
        &data3[15*n],
        &data4[15*n],
        &data5[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[49].time     += wall_t2 - wall_t1;
  OP_kernels[49].transfer += (float)set->size * arg1.size;
  OP_kernels[49].transfer += (float)set->size * arg2.size;
  OP_kernels[49].transfer += (float)set->size * arg3.size;
  OP_kernels[49].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[49].transfer += (float)set->size * arg5.size * 2.0f;
}
