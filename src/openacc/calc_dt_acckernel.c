//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void calc_dt_openacc( const double *x, const double *y, double *dt) {
  double len0 = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]));
  double len1 = sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]));
  double len2 = sqrt((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0]));
  double sper = (len0 + len1 + len2) / 2.0;
  double area = sqrt(sper * (sper - len0) * (sper - len1) * (sper - len2));
  if(*dt > area / sper)
    *dt = area / sper;
}

// host stub function
void op_par_loop_calc_dt(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  double*arg2h = (double *)arg2.data;
  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(45);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[45].name      = name;
  OP_kernels[45].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  calc_dt");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg2_l = arg2h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1) reduction(min:arg2_l)
    for ( int n=0; n<set->size; n++ ){
      calc_dt_openacc(
        &data0[3*n],
        &data1[3*n],
        &arg2_l);
    }
  }

  // combine reduction data
  arg2h[0]  = MIN(arg2h[0],arg2_l);
  op_mpi_reduce_double(&arg2,arg2h);
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[45].time     += wall_t2 - wall_t1;
  OP_kernels[45].transfer += (float)set->size * arg0.size;
  OP_kernels[45].transfer += (float)set->size * arg1.size;
}
