//
// auto-generated by op2.py
//

//user function
//user function

void ls_reinit_check_omp4_kernel(
  double *arg0,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  double *arg4,
  int *arg5,
  int count,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_ls_reinit_check(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  double*arg0h = (double *)arg0.data;
  double*arg4h = (double *)arg4.data;
  int*arg5h = (int *)arg5.data;
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
  op_timing_realloc(74);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[74].name      = name;
  OP_kernels[74].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_reinit_check");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_74
    int part_size = OP_PART_SIZE_74;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_74
    int nthread = OP_BLOCK_SIZE_74;
  #else
    int nthread = OP_block_size;
  #endif

  double arg0_l = arg0h[0];
  double arg4_l = arg4h[0];
  int arg5_l = arg5h[0];

  if (set_size >0) {

    //Set up typed device pointers for OpenMP

    double* data1 = (double*)arg1.data_d;
    int dat1size = getSetSizeFromOpArg(&arg1) * arg1.dat->dim;
    double* data2 = (double*)arg2.data_d;
    int dat2size = getSetSizeFromOpArg(&arg2) * arg2.dat->dim;
    double* data3 = (double*)arg3.data_d;
    int dat3size = getSetSizeFromOpArg(&arg3) * arg3.dat->dim;
    ls_reinit_check_omp4_kernel(
      &arg0_l,
      data1,
      dat1size,
      data2,
      dat2size,
      data3,
      dat3size,
      &arg4_l,
      &arg5_l,
      set->size,
      part_size!=0?(set->size-1)/part_size+1:(set->size-1)/nthread,
      nthread);

  }

  // combine reduction data
  arg4h[0] = arg4_l;
  op_mpi_reduce_double(&arg4,arg4h);
  arg5h[0] = arg5_l;
  op_mpi_reduce_int(&arg5,arg5h);
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[74].time     += wall_t2 - wall_t1;
  OP_kernels[74].transfer += (float)set->size * arg1.size;
  OP_kernels[74].transfer += (float)set->size * arg2.size;
  OP_kernels[74].transfer += (float)set->size * arg3.size;
}
