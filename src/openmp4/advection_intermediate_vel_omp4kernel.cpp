//
// auto-generated by op2.py
//

//user function
//user function

void advection_intermediate_vel_omp4_kernel(
  double *arg0,
  double *arg1,
  double *arg2,
  double *arg3,
  double *arg4,
  double *arg5,
  double *data6,
  int dat6size,
  double *data7,
  int dat7size,
  double *data8,
  int dat8size,
  double *data9,
  int dat9size,
  double *data10,
  int dat10size,
  double *data11,
  int dat11size,
  double *data12,
  int dat12size,
  double *data13,
  int dat13size,
  double *data14,
  int dat14size,
  double *data15,
  int dat15size,
  int count,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_advection_intermediate_vel(char const *name, op_set set,
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
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15){

  double*arg0h = (double *)arg0.data;
  double*arg1h = (double *)arg1.data;
  double*arg2h = (double *)arg2.data;
  double*arg3h = (double *)arg3.data;
  double*arg4h = (double *)arg4.data;
  double*arg5h = (double *)arg5.data;
  int nargs = 16;
  op_arg args[16];

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
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(6);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[6].name      = name;
  OP_kernels[6].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  advection_intermediate_vel");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_6
    int part_size = OP_PART_SIZE_6;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_6
    int nthread = OP_BLOCK_SIZE_6;
  #else
    int nthread = OP_block_size;
  #endif

  double arg0_l = arg0h[0];
  double arg1_l = arg1h[0];
  double arg2_l = arg2h[0];
  double arg3_l = arg3h[0];
  double arg4_l = arg4h[0];
  double arg5_l = arg5h[0];

  if (set_size >0) {

    //Set up typed device pointers for OpenMP

    double* data6 = (double*)arg6.data_d;
    int dat6size = getSetSizeFromOpArg(&arg6) * arg6.dat->dim;
    double* data7 = (double*)arg7.data_d;
    int dat7size = getSetSizeFromOpArg(&arg7) * arg7.dat->dim;
    double* data8 = (double*)arg8.data_d;
    int dat8size = getSetSizeFromOpArg(&arg8) * arg8.dat->dim;
    double* data9 = (double*)arg9.data_d;
    int dat9size = getSetSizeFromOpArg(&arg9) * arg9.dat->dim;
    double* data10 = (double*)arg10.data_d;
    int dat10size = getSetSizeFromOpArg(&arg10) * arg10.dat->dim;
    double* data11 = (double*)arg11.data_d;
    int dat11size = getSetSizeFromOpArg(&arg11) * arg11.dat->dim;
    double* data12 = (double*)arg12.data_d;
    int dat12size = getSetSizeFromOpArg(&arg12) * arg12.dat->dim;
    double* data13 = (double*)arg13.data_d;
    int dat13size = getSetSizeFromOpArg(&arg13) * arg13.dat->dim;
    double* data14 = (double*)arg14.data_d;
    int dat14size = getSetSizeFromOpArg(&arg14) * arg14.dat->dim;
    double* data15 = (double*)arg15.data_d;
    int dat15size = getSetSizeFromOpArg(&arg15) * arg15.dat->dim;
    advection_intermediate_vel_omp4_kernel(
      &arg0_l,
      &arg1_l,
      &arg2_l,
      &arg3_l,
      &arg4_l,
      &arg5_l,
      data6,
      dat6size,
      data7,
      dat7size,
      data8,
      dat8size,
      data9,
      dat9size,
      data10,
      dat10size,
      data11,
      dat11size,
      data12,
      dat12size,
      data13,
      dat13size,
      data14,
      dat14size,
      data15,
      dat15size,
      set->size,
      part_size!=0?(set->size-1)/part_size+1:(set->size-1)/nthread,
      nthread);

  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[6].time     += wall_t2 - wall_t1;
  OP_kernels[6].transfer += (float)set->size * arg6.size;
  OP_kernels[6].transfer += (float)set->size * arg7.size;
  OP_kernels[6].transfer += (float)set->size * arg8.size;
  OP_kernels[6].transfer += (float)set->size * arg9.size;
  OP_kernels[6].transfer += (float)set->size * arg10.size;
  OP_kernels[6].transfer += (float)set->size * arg11.size;
  OP_kernels[6].transfer += (float)set->size * arg12.size;
  OP_kernels[6].transfer += (float)set->size * arg13.size;
  OP_kernels[6].transfer += (float)set->size * arg14.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg15.size * 2.0f;
}
