//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void set_rkQ_openacc( const int *stage, const double *dt, const double *s,
                    const double *rk0, const double *rk1, double *rkQ) {
  if(*stage == -1) {
    for(int i = 0; i < 10; i++) {
      rkQ[i] = s[i];
    }
  } else if(*stage == 0) {
    for(int i = 0; i < 10; i++) {
      rkQ[i] = s[i] + (*dt) * rk0[i];
    }
  } else {
    for(int i = 0; i < 10; i++) {
      rkQ[i] = s[i] + 0.5 * (*dt) * (rk0[i] / 4.0 + rk1[i] / 4.0);
    }
  }
}

// host stub function
void op_par_loop_set_rkQ(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int*arg0h = (int *)arg0.data;
  double*arg1h = (double *)arg1.data;
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
  op_timing_realloc(47);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[47].name      = name;
  OP_kernels[47].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  set_rkQ");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  int arg0_l = arg0h[0];
  double arg1_l = arg1h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    #pragma acc parallel loop independent deviceptr(data2,data3,data4,data5)
    for ( int n=0; n<set->size; n++ ){
      set_rkQ_openacc(
        &arg0_l,
        &arg1_l,
        &data2[10*n],
        &data3[10*n],
        &data4[10*n],
        &data5[10*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[47].time     += wall_t2 - wall_t1;
  OP_kernels[47].transfer += (float)set->size * arg2.size;
  OP_kernels[47].transfer += (float)set->size * arg3.size;
  OP_kernels[47].transfer += (float)set->size * arg4.size;
  OP_kernels[47].transfer += (float)set->size * arg5.size * 2.0f;
}
