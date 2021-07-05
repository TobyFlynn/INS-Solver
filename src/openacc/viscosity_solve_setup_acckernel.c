//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void viscosity_solve_setup_openacc( const double *mu, const double *rho,
                                  const double *mmConst, double *factor,
                                  double *mmFactor) {
  for(int i = 0; i < 15; i++) {
    factor[i] = mu[i];
    mmFactor[i] = *mmConst * rho[i];
  }
}

// host stub function
void op_par_loop_viscosity_solve_setup(char const *name, op_set set,
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
  op_timing_realloc(24);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[24].name      = name;
  OP_kernels[24].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  viscosity_solve_setup");
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
      viscosity_solve_setup_openacc(
        &data0[15*n],
        &data1[15*n],
        &arg2_l,
        &data3[15*n],
        &data4[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[24].time     += wall_t2 - wall_t1;
  OP_kernels[24].transfer += (float)set->size * arg0.size;
  OP_kernels[24].transfer += (float)set->size * arg1.size;
  OP_kernels[24].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[24].transfer += (float)set->size * arg4.size * 2.0f;
}
