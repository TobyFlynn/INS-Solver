//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void ls_step_openacc( const double *alpha, const double *s, double *step,
                    double *nu) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 15; i++) {
    step[i] = tanh(PI * s[i] / *alpha);
    nu[i] = nu0 * step[i] + nu1 * (1.0 - step[i]);
  }
}

// host stub function
void op_par_loop_ls_step(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  double*arg0h = (double *)arg0.data;
  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(73);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[73].name      = name;
  OP_kernels[73].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_step");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg0_l = arg0h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    #pragma acc parallel loop independent deviceptr(data1,data2,data3)
    for ( int n=0; n<set->size; n++ ){
      ls_step_openacc(
        &arg0_l,
        &data1[15*n],
        &data2[15*n],
        &data3[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[73].time     += wall_t2 - wall_t1;
  OP_kernels[73].transfer += (float)set->size * arg1.size;
  OP_kernels[73].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[73].transfer += (float)set->size * arg3.size * 2.0f;
}
