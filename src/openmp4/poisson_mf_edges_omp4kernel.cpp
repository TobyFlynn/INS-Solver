//
// auto-generated by op2.py
//

//user function
//user function

void poisson_mf_edges_omp4_kernel(
  int *data0,
  int dat0size,
  bool *data1,
  int dat1size,
  int *map2,
  int map2size,
  double *data2,
  int dat2size,
  double *data4,
  int dat4size,
  double *data6,
  int dat6size,
  double *data8,
  int dat8size,
  double *data10,
  int dat10size,
  double *data12,
  int dat12size,
  double *data14,
  int dat14size,
  double *data16,
  int dat16size,
  double *data18,
  int dat18size,
  double *data20,
  int dat20size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_poisson_mf_edges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6,
  op_arg arg8,
  op_arg arg10,
  op_arg arg12,
  op_arg arg14,
  op_arg arg16,
  op_arg arg18,
  op_arg arg20){

  int nargs = 22;
  op_arg args[22];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 21, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 21, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 21, "double", OP_READ);
  }

  arg8.idx = 0;
  args[8] = arg8;
  for ( int v=1; v<2; v++ ){
    args[8 + v] = op_arg_dat(arg8.dat, v, arg8.map, 3, "double", OP_READ);
  }

  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 21, "double", OP_READ);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 21, "double", OP_READ);
  }

  arg14.idx = 0;
  args[14] = arg14;
  for ( int v=1; v<2; v++ ){
    args[14 + v] = op_arg_dat(arg14.dat, v, arg14.map, 21, "double", OP_READ);
  }

  arg16.idx = 0;
  args[16] = arg16;
  for ( int v=1; v<2; v++ ){
    args[16 + v] = op_arg_dat(arg16.dat, v, arg16.map, 21, "double", OP_INC);
  }

  arg18.idx = 0;
  args[18] = arg18;
  for ( int v=1; v<2; v++ ){
    args[18 + v] = op_arg_dat(arg18.dat, v, arg18.map, 21, "double", OP_INC);
  }

  arg20.idx = 0;
  args[20] = arg20;
  for ( int v=1; v<2; v++ ){
    args[20 + v] = op_arg_dat(arg20.dat, v, arg20.map, 21, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(31);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[31].name      = name;
  OP_kernels[31].count    += 1;

  int  ninds   = 10;
  int  inds[22] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_mf_edges\n");
  }

  // get plan
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_31
    int part_size = OP_PART_SIZE_31;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_31
    int nthread = OP_BLOCK_SIZE_31;
  #else
    int nthread = OP_block_size;
  #endif


  int ncolors = 0;
  int set_size1 = set->size + set->exec_size;

  if (set_size >0) {

    //Set up typed device pointers for OpenMP
    int *map2 = arg2.map_data_d;
     int map2size = arg2.map->dim * set_size1;

    int* data0 = (int*)arg0.data_d;
    int dat0size = getSetSizeFromOpArg(&arg0) * arg0.dat->dim;
    bool* data1 = (bool*)arg1.data_d;
    int dat1size = getSetSizeFromOpArg(&arg1) * arg1.dat->dim;
    double *data2 = (double *)arg2.data_d;
    int dat2size = getSetSizeFromOpArg(&arg2) * arg2.dat->dim;
    double *data4 = (double *)arg4.data_d;
    int dat4size = getSetSizeFromOpArg(&arg4) * arg4.dat->dim;
    double *data6 = (double *)arg6.data_d;
    int dat6size = getSetSizeFromOpArg(&arg6) * arg6.dat->dim;
    double *data8 = (double *)arg8.data_d;
    int dat8size = getSetSizeFromOpArg(&arg8) * arg8.dat->dim;
    double *data10 = (double *)arg10.data_d;
    int dat10size = getSetSizeFromOpArg(&arg10) * arg10.dat->dim;
    double *data12 = (double *)arg12.data_d;
    int dat12size = getSetSizeFromOpArg(&arg12) * arg12.dat->dim;
    double *data14 = (double *)arg14.data_d;
    int dat14size = getSetSizeFromOpArg(&arg14) * arg14.dat->dim;
    double *data16 = (double *)arg16.data_d;
    int dat16size = getSetSizeFromOpArg(&arg16) * arg16.dat->dim;
    double *data18 = (double *)arg18.data_d;
    int dat18size = getSetSizeFromOpArg(&arg18) * arg18.dat->dim;
    double *data20 = (double *)arg20.data_d;
    int dat20size = getSetSizeFromOpArg(&arg20) * arg20.dat->dim;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      poisson_mf_edges_omp4_kernel(
        data0,
        dat0size,
        data1,
        dat1size,
        map2,
        map2size,
        data2,
        dat2size,
        data4,
        dat4size,
        data6,
        dat6size,
        data8,
        dat8size,
        data10,
        dat10size,
        data12,
        dat12size,
        data14,
        dat14size,
        data16,
        dat16size,
        data18,
        dat18size,
        data20,
        dat20size,
        col_reord,
        set_size1,
        start,
        end,
        part_size!=0?(end-start-1)/part_size+1:(end-start-1)/nthread,
        nthread);

    }
    OP_kernels[31].transfer  += Plan->transfer;
    OP_kernels[31].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[31].time     += wall_t2 - wall_t1;
}
