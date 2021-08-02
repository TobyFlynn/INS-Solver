//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void ls_reinit_check_openacc( const double *alpha, const double *s,
                            const double *dsdx, const double *dsdy,
                            double *res, int *count) {
  for(int i = 0; i < 10; i++) {
    if(fabs(s[i]) < (*alpha)) {
      *res += dsdx[i] * dsdx[i] + dsdy[i] * dsdy[i];
      *count += 1;
    }
  }
}

// host stub function
void op_par_loop_ls_reinit_check(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  double*arg0h = (double *)arg0.data;
  double*arg4h = (double *)arg4.data;
  int*arg5h = (int *)arg5.data;
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
  op_timing_realloc(63);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[63].name      = name;
  OP_kernels[63].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_reinit_check");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg0_l = arg0h[0];
  double arg4_l = arg4h[0];
  int arg5_l = arg5h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    #pragma acc parallel loop independent deviceptr(data1,data2,data3) reduction(+:arg4_l) reduction(+:arg5_l)
    for ( int n=0; n<set->size; n++ ){
      ls_reinit_check_openacc(
        &arg0_l,
        &data1[10*n],
        &data2[10*n],
        &data3[10*n],
        &arg4_l,
        &arg5_l);
    }
  }

  // combine reduction data
  arg4h[0] = arg4_l;
  op_mpi_reduce_double(&arg4,arg4h);
  arg5h[0] = arg5_l;
  op_mpi_reduce_int(&arg5,arg5h);
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[63].time     += wall_t2 - wall_t1;
  OP_kernels[63].transfer += (float)set->size * arg1.size;
  OP_kernels[63].transfer += (float)set->size * arg2.size;
  OP_kernels[63].transfer += (float)set->size * arg3.size;
}
