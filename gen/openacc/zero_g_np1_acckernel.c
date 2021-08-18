//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void zero_g_np1_openacc( double *g_np0) {
  for(int i = 0; i < 18; i++) {
    g_np0[i] = 0.0;
  }
}

// host stub function
void op_par_loop_zero_g_np1(char const *name, op_set set,
  op_arg arg0){

  int nargs = 1;
  op_arg args[1];

  args[0] = arg0;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(35);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[35].name      = name;
  OP_kernels[35].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  zero_g_np1");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    #pragma acc parallel loop independent deviceptr(data0)
    for ( int n=0; n<set->size; n++ ){
      zero_g_np1_openacc(
        &data0[18*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[35].time     += wall_t2 - wall_t1;
  OP_kernels[35].transfer += (float)set->size * arg0.size * 2.0f;
}
