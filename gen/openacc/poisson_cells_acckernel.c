//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_cells_openacc( const double *u, const double *op, double *rhs) {
  for(int m = 0; m < 3; m++) {
    int ind = m * 3;
    rhs[m] = 0.0;
    for(int n = 0; n < 3; n++) {
      rhs[m] += op[ind + n] * u[n];
    }
  }
}

// host stub function
void op_par_loop_poisson_cells(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(16);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[16].name      = name;
  OP_kernels[16].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_cells");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2)
    for ( int n=0; n<set->size; n++ ){
      poisson_cells_openacc(
        &data0[3*n],
        &data1[9*n],
        &data2[3*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[16].time     += wall_t2 - wall_t1;
  OP_kernels[16].transfer += (float)set->size * arg0.size;
  OP_kernels[16].transfer += (float)set->size * arg1.size;
  OP_kernels[16].transfer += (float)set->size * arg2.size * 2.0f;
}
