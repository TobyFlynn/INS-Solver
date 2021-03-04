//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void init_cubature_OP_openacc( const double *J, const double *Dx, const double *Dy,
                             double *temp, double *temp2) {
  for(int m = 0; m < 46; m++) {
    for(int n = 0; n < 15; n++) {
      int ind = m * 15 + n;
      temp[ind] = J[m] * cubW[m] * Dx[ind];
      temp2[ind] = J[m] * cubW[m] * Dy[ind];
    }
  }
}

// host stub function
void op_par_loop_init_cubature_OP(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(15);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[15].name      = name;
  OP_kernels[15].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  init_cubature_OP");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3,data4)
    for ( int n=0; n<set->size; n++ ){
      init_cubature_OP_openacc(
        &data0[46*n],
        &data1[690*n],
        &data2[690*n],
        &data3[690*n],
        &data4[690*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[15].time     += wall_t2 - wall_t1;
  OP_kernels[15].transfer += (float)set->size * arg0.size;
  OP_kernels[15].transfer += (float)set->size * arg1.size;
  OP_kernels[15].transfer += (float)set->size * arg2.size;
  OP_kernels[15].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[15].transfer += (float)set->size * arg4.size * 2.0f;
}
