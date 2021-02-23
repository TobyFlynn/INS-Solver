//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void viscosity_rhs_openacc( const double *factor, const double *J, const double *qtt0,
                          const double *qtt1, double *vRHS0, double *vRHS1) {

  for(int i = 0; i < 15; i++) {
    vRHS0[i] = *factor * J[i] * qtt0[i];
    vRHS1[i] = *factor * J[i] * qtt1[i];


  }
}

// host stub function
void op_par_loop_viscosity_rhs(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  double*arg0h = (double *)arg0.data;
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
  op_timing_realloc(15);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[15].name      = name;
  OP_kernels[15].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  viscosity_rhs");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  double arg0_l = arg0h[0];

  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    #pragma acc parallel loop independent deviceptr(data1,data2,data3,data4,data5)
    for ( int n=0; n<set->size; n++ ){
      viscosity_rhs_openacc(
        &arg0_l,
        &data1[15*n],
        &data2[15*n],
        &data3[15*n],
        &data4[15*n],
        &data5[15*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[15].time     += wall_t2 - wall_t1;
  OP_kernels[15].transfer += (float)set->size * arg1.size;
  OP_kernels[15].transfer += (float)set->size * arg2.size;
  OP_kernels[15].transfer += (float)set->size * arg3.size;
  OP_kernels[15].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[15].transfer += (float)set->size * arg5.size * 2.0f;
}
