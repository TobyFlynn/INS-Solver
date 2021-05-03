//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_mf2_mass_openacc( const double *u, const double *op, const double *factor,
                             const double *mm, double *rhs) {
  double mFactor = *factor;
  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    double val = 0.0;
    for(int n = 0; n < 15; n++) {

      val += (op[ind + n] + mm[n * 15 + m] * mFactor) * u[n];
    }
    rhs[m] = val;
  }
}

// host stub function
void op_par_loop_poisson_mf2_mass(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  double*arg2h = (double *)arg2.data;
  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(49);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[49].name      = name;
  OP_kernels[49].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_mf2_mass");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg2_l = arg2h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data3,data4)
    for ( int n=0; n<set->size; n++ ){
      poisson_mf2_mass_openacc(
        &data0[15*n],
        &data1[225*n],
        &arg2_l,
        &data3[225*n],
        &data4[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[49].time     += wall_t2 - wall_t1;
  OP_kernels[49].transfer += (float)set->size * arg0.size;
  OP_kernels[49].transfer += (float)set->size * arg1.size;
  OP_kernels[49].transfer += (float)set->size * arg3.size;
  OP_kernels[49].transfer += (float)set->size * arg4.size * 2.0f;
}