//
// auto-generated by op2.py
//

//user function
//user function

void pressure_rhs_omp4_kernel(
  double *arg0,
  double *arg1,
  double *arg2,
  double *arg3,
  double *data4,
  int dat4size,
  double *data5,
  int dat5size,
  double *data6,
  int dat6size,
  double *data7,
  int dat7size,
  double *data8,
  int dat8size,
  int count,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_pressure_rhs(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  double*arg0h = (double *)arg0.data;
  double*arg1h = (double *)arg1.data;
  double*arg2h = (double *)arg2.data;
  double*arg3h = (double *)arg3.data;
  int nargs = 9;
  op_arg args[9];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(58);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[58].name      = name;
  OP_kernels[58].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pressure_rhs");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_58
    int part_size = OP_PART_SIZE_58;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_58
    int nthread = OP_BLOCK_SIZE_58;
  #else
    int nthread = OP_block_size;
  #endif

  double arg0_l = arg0h[0];
  double arg1_l = arg1h[0];
  double arg2_l = arg2h[0];
  double arg3_l = arg3h[0];

  if (set_size >0) {

    //Set up typed device pointers for OpenMP

    double* data4 = (double*)arg4.data_d;
    int dat4size = getSetSizeFromOpArg(&arg4) * arg4.dat->dim;
    double* data5 = (double*)arg5.data_d;
    int dat5size = getSetSizeFromOpArg(&arg5) * arg5.dat->dim;
    double* data6 = (double*)arg6.data_d;
    int dat6size = getSetSizeFromOpArg(&arg6) * arg6.dat->dim;
    double* data7 = (double*)arg7.data_d;
    int dat7size = getSetSizeFromOpArg(&arg7) * arg7.dat->dim;
    double* data8 = (double*)arg8.data_d;
    int dat8size = getSetSizeFromOpArg(&arg8) * arg8.dat->dim;
    pressure_rhs_omp4_kernel(
      &arg0_l,
      &arg1_l,
      &arg2_l,
      &arg3_l,
      data4,
      dat4size,
      data5,
      dat5size,
      data6,
      dat6size,
      data7,
      dat7size,
      data8,
      dat8size,
      set->size,
      part_size!=0?(set->size-1)/part_size+1:(set->size-1)/nthread,
      nthread);

  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[58].time     += wall_t2 - wall_t1;
  OP_kernels[58].transfer += (float)set->size * arg4.size;
  OP_kernels[58].transfer += (float)set->size * arg5.size;
  OP_kernels[58].transfer += (float)set->size * arg6.size;
  OP_kernels[58].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[58].transfer += (float)set->size * arg8.size * 2.0f;
}
