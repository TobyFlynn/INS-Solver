//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_bc2_openacc( const double *sJ, const double *tau, const double *bc,
                        double *gtau) {
  for(int i = 0; i < 21; i++) {
    gtau[i] = tau[i / 7] * gaussW_g[i % 7] * sJ[i] * bc[i];
  }
}

// host stub function
void op_par_loop_poisson_bc2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(46);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[46].name      = name;
  OP_kernels[46].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_bc2");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3)
    for ( int n=0; n<set->size; n++ ){
      poisson_bc2_openacc(
        &data0[21*n],
        &data1[3*n],
        &data2[21*n],
        &data3[21*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[46].time     += wall_t2 - wall_t1;
  OP_kernels[46].transfer += (float)set->size * arg0.size;
  OP_kernels[46].transfer += (float)set->size * arg1.size;
  OP_kernels[46].transfer += (float)set->size * arg2.size;
  OP_kernels[46].transfer += (float)set->size * arg3.size * 2.0f;
}
