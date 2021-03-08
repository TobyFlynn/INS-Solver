//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_test_error_openacc( const double *x, const double *y,
                               const double *sol, double *err, double *l2) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 15; i++) {
    double x1 = x[i];
    double y1 = y[i];


    double exact = (1.0 - (x[i] * x[i])) * (2.0 * (y[i] * y[i] * y[i]) - 3.0 * (y[i] * y[i]) + 1.0);

    err[i] = fabs(sol[i] - exact);
    *l2 += err[i] * err[i];
  }
}

// host stub function
void op_par_loop_poisson_test_error(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  double*arg4h = (double *)arg4.data;
  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(33);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[33].name      = name;
  OP_kernels[33].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_test_error");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg4_l = arg4h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3) reduction(+:arg4_l)
    for ( int n=0; n<set->size; n++ ){
      poisson_test_error_openacc(
        &data0[15*n],
        &data1[15*n],
        &data2[15*n],
        &data3[15*n],
        &arg4_l);
    }
  }

  // combine reduction data
  arg4h[0] = arg4_l;
  op_mpi_reduce_double(&arg4,arg4h);
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[33].time     += wall_t2 - wall_t1;
  OP_kernels[33].transfer += (float)set->size * arg0.size;
  OP_kernels[33].transfer += (float)set->size * arg1.size;
  OP_kernels[33].transfer += (float)set->size * arg2.size;
  OP_kernels[33].transfer += (float)set->size * arg3.size * 2.0f;
}
