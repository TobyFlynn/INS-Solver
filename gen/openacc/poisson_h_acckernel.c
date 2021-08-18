//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_h_openacc( const double *x, const double *y, double *h) {
  double len0 = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]));
  double len1 = sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]));
  double len2 = sqrt((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0]));
  double sper = (len0 + len1 + len2) / 2.0;
  double area = sqrt(sper * (sper - len0) * (sper - len1) * (sper - len2));
  *h = 2.0 * sper / area;
}

// host stub function
void op_par_loop_poisson_h(char const *name, op_set set,
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
  op_timing_realloc(17);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[17].name      = name;
  OP_kernels[17].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_h");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2)
    for ( int n=0; n<set->size; n++ ){
      poisson_h_openacc(
        &data0[3*n],
        &data1[3*n],
        &data2[1*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[17].time     += wall_t2 - wall_t1;
  OP_kernels[17].transfer += (float)set->size * arg0.size;
  OP_kernels[17].transfer += (float)set->size * arg1.size;
  OP_kernels[17].transfer += (float)set->size * arg2.size * 2.0f;
}
