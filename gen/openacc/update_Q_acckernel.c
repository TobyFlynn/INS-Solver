//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void update_Q_openacc( const double *dt, double *s, const double *rk0,
                     const double *rk1, const double *rk2) {
  for(int i = 0; i < 10; i++) {
    double add = (*dt) * (rk0[i]/ 6.0 + rk1[i] / 6.0 + 2.0 * rk2[i] / 3.0);
    s[i] = s[i] + add;
  }
}

// host stub function
void op_par_loop_update_Q(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  double*arg0h = (double *)arg0.data;
  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(12);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[12].name      = name;
  OP_kernels[12].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  update_Q");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg0_l = arg0h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    #pragma acc parallel loop independent deviceptr(data1,data2,data3,data4)
    for ( int n=0; n<set->size; n++ ){
      update_Q_openacc(
        &arg0_l,
        &data1[10*n],
        &data2[10*n],
        &data3[10*n],
        &data4[10*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[12].time     += wall_t2 - wall_t1;
  OP_kernels[12].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[12].transfer += (float)set->size * arg2.size;
  OP_kernels[12].transfer += (float)set->size * arg3.size;
  OP_kernels[12].transfer += (float)set->size * arg4.size;
}
