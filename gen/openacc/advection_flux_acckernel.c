//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void advection_flux_openacc( const double *u, const double *v, double *f0,
                           double *f1, double *f2, double *f3) {
  for(int i = 0; i < 3; i++) {
    f0[i] = u[i] * u[i];
    f1[i] = u[i] * v[i];
    f2[i] = u[i] * v[i];
    f3[i] = v[i] * v[i];
  }
}

// host stub function
void op_par_loop_advection_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

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
  op_timing_realloc(28);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[28].name      = name;
  OP_kernels[28].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  advection_flux");
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
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4,data5)
    for ( int n=0; n<set->size; n++ ){
      advection_flux_openacc(
        &data0[3*n],
        &data1[3*n],
        &data2[3*n],
        &data3[3*n],
        &data4[3*n],
        &data5[3*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[28].time     += wall_t2 - wall_t1;
  OP_kernels[28].transfer += (float)set->size * arg0.size;
  OP_kernels[28].transfer += (float)set->size * arg1.size;
  OP_kernels[28].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[28].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[28].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[28].transfer += (float)set->size * arg5.size * 2.0f;
}
