//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_mf2_mass_openacc( const double *u, const double *factor,
                             const double *mm, double *rhs) {
  double mFactor = *factor;
  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {

      val += (mm[n * 15 + m] * mFactor) * u[n];
    }
    rhs[m] = val;
  }
}

// host stub function
void op_par_loop_poisson_mf2_mass(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  double*arg1h = (double *)arg1.data;
  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(38);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[38].name      = name;
  OP_kernels[38].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_mf2_mass");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg1_l = arg1h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data2,data3)
    for ( int n=0; n<set->size; n++ ){
      poisson_mf2_mass_openacc(
        &data0[15*n],
        &arg1_l,
        &data2[225*n],
        &data3[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[38].time     += wall_t2 - wall_t1;
  OP_kernels[38].transfer += (float)set->size * arg0.size;
  OP_kernels[38].transfer += (float)set->size * arg2.size;
  OP_kernels[38].transfer += (float)set->size * arg3.size * 2.0f;
}
