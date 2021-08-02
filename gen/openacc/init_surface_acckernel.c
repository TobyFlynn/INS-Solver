//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void init_surface_openacc( const double *x, const double *y, double *s) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 10; i++) {

    s[i] = sqrt((x[i] - 1.0) * (x[i] - 1.0) + (y[i] - 0.5) * (y[i] - 0.5)) - 0.1;
  }
}

// host stub function
void op_par_loop_init_surface(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(46);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[46].name      = name;
  OP_kernels[46].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_surface");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2)
    for ( int n=0; n<set->size; n++ ){
      init_surface_openacc(
        &data0[10*n],
        &data1[10*n],
        &data2[10*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[46].time     += wall_t2 - wall_t1;
  OP_kernels[46].transfer += (float)set->size * arg0.size;
  OP_kernels[46].transfer += (float)set->size * arg1.size;
  OP_kernels[46].transfer += (float)set->size * arg2.size * 2.0f;
}
