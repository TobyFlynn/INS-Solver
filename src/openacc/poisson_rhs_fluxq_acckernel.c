//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_rhs_fluxq_openacc( const double *nx, const double *ny, const double *fscale,
                       const double *tau, const double *u, const double *du, const double *qx,
                       const double *qy, double *exQx, double *exQy, double *fluxq) {
  for(int i = 0; i < 15; i++) {
    double dqx = (qx[FMASK[i]] + exQx[i]) / 2.0;
    double dqy = (qy[FMASK[i]] + exQy[i]) / 2.0;
    fluxq[i] = fscale[i] * (nx[i] * dqx + ny[i] * dqy - tau[i] * (u[FMASK[i]] - du[i]));
    exQx[i] = 0.0;
    exQy[i] = 0.0;
  }
}

// host stub function
void op_par_loop_poisson_rhs_fluxq(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10){

  int nargs = 11;
  op_arg args[11];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  args[10] = arg10;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(27);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[27].name      = name;
  OP_kernels[27].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_rhs_fluxq");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    double* data6 = (double*)arg6.data_d;
    double* data7 = (double*)arg7.data_d;
    double* data8 = (double*)arg8.data_d;
    double* data9 = (double*)arg9.data_d;
    double* data10 = (double*)arg10.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4,data5,data6,data7,data8,data9,data10)
    for ( int n=0; n<set->size; n++ ){
      poisson_rhs_fluxq_openacc(
        &data0[15*n],
        &data1[15*n],
        &data2[15*n],
        &data3[15*n],
        &data4[15*n],
        &data5[15*n],
        &data6[15*n],
        &data7[15*n],
        &data8[15*n],
        &data9[15*n],
        &data10[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[27].time     += wall_t2 - wall_t1;
  OP_kernels[27].transfer += (float)set->size * arg0.size;
  OP_kernels[27].transfer += (float)set->size * arg1.size;
  OP_kernels[27].transfer += (float)set->size * arg2.size;
  OP_kernels[27].transfer += (float)set->size * arg3.size;
  OP_kernels[27].transfer += (float)set->size * arg4.size;
  OP_kernels[27].transfer += (float)set->size * arg5.size;
  OP_kernels[27].transfer += (float)set->size * arg6.size;
  OP_kernels[27].transfer += (float)set->size * arg7.size;
  OP_kernels[27].transfer += (float)set->size * arg8.size * 2.0f;
  OP_kernels[27].transfer += (float)set->size * arg9.size * 2.0f;
  OP_kernels[27].transfer += (float)set->size * arg10.size * 2.0f;
}
