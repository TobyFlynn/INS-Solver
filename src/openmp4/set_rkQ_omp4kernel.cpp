//
// auto-generated by op2.py
//

//user function
//user function

void set_rkQ_omp4_kernel(
  int *arg0,
  double *arg1,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  double *data5,
  int dat5size,
  int count,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_set_rkQ(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int*arg0h = (int *)arg0.data;
  double*arg1h = (double *)arg1.data;
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
  op_timing_realloc(54);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[54].name      = name;
  OP_kernels[54].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  set_rkQ");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_54
    int part_size = OP_PART_SIZE_54;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_54
    int nthread = OP_BLOCK_SIZE_54;
  #else
    int nthread = OP_block_size;
  #endif

  int arg0_l = arg0h[0];
  double arg1_l = arg1h[0];

  if (set_size >0) {

    //Set up typed device pointers for OpenMP

    double* data2 = (double*)arg2.data_d;
    int dat2size = getSetSizeFromOpArg(&arg2) * arg2.dat->dim;
    double* data3 = (double*)arg3.data_d;
    int dat3size = getSetSizeFromOpArg(&arg3) * arg3.dat->dim;
    double* data4 = (double*)arg4.data_d;
    int dat4size = getSetSizeFromOpArg(&arg4) * arg4.dat->dim;
    double* data5 = (double*)arg5.data_d;
    int dat5size = getSetSizeFromOpArg(&arg5) * arg5.dat->dim;
    set_rkQ_omp4_kernel(
      &arg0_l,
      &arg1_l,
      data2,
      dat2size,
      data3,
      dat3size,
      data4,
      dat4size,
      data5,
      dat5size,
      set->size,
      part_size!=0?(set->size-1)/part_size+1:(set->size-1)/nthread,
      nthread);

  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[54].time     += wall_t2 - wall_t1;
  OP_kernels[54].transfer += (float)set->size * arg2.size;
  OP_kernels[54].transfer += (float)set->size * arg3.size;
  OP_kernels[54].transfer += (float)set->size * arg4.size;
  OP_kernels[54].transfer += (float)set->size * arg5.size * 2.0f;
}
