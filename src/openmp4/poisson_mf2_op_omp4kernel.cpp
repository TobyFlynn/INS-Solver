//
// auto-generated by op2.py
//

//user function
//user function

void poisson_mf2_op_omp4_kernel(
  double *data0,
  int dat0size,
  double *arg1,
  double *data2,
  int dat2size,
  int count,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_poisson_mf2_op(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  double*arg1h = (double *)arg1.data;
  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(21);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[21].name      = name;
  OP_kernels[21].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_mf2_op");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_21
    int part_size = OP_PART_SIZE_21;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_21
    int nthread = OP_BLOCK_SIZE_21;
  #else
    int nthread = OP_block_size;
  #endif

  double arg1_l = arg1h[0];

  if (set_size >0) {

    //Set up typed device pointers for OpenMP

    double* data0 = (double*)arg0.data_d;
    int dat0size = getSetSizeFromOpArg(&arg0) * arg0.dat->dim;
    double* data2 = (double*)arg2.data_d;
    int dat2size = getSetSizeFromOpArg(&arg2) * arg2.dat->dim;
    poisson_mf2_op_omp4_kernel(
      data0,
      dat0size,
      &arg1_l,
      data2,
      dat2size,
      set->size,
      part_size!=0?(set->size-1)/part_size+1:(set->size-1)/nthread,
      nthread);

  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[21].time     += wall_t2 - wall_t1;
  OP_kernels[21].transfer += (float)set->size * arg0.size;
  OP_kernels[21].transfer += (float)set->size * arg2.size * 2.0f;
}
