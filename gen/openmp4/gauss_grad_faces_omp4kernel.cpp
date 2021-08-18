//
// auto-generated by op2.py
//

//user function
//user function

void gauss_grad_faces_omp4_kernel(
  int *data0,
  int dat0size,
  int *map1,
  int map1size,
  double *data1,
  int dat1size,
  double *data3,
  int dat3size,
  double *data5,
  int dat5size,
  double *data7,
  int dat7size,
  double *data9,
  int dat9size,
  double *data11,
  int dat11size,
  double *data13,
  int dat13size,
  double *data15,
  int dat15size,
  double *data17,
  int dat17size,
  double *data19,
  int dat19size,
  double *data21,
  int dat21size,
  double *data23,
  int dat23size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_gauss_grad_faces(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3,
  op_arg arg5,
  op_arg arg7,
  op_arg arg9,
  op_arg arg11,
  op_arg arg13,
  op_arg arg15,
  op_arg arg17,
  op_arg arg19,
  op_arg arg21,
  op_arg arg23){

  int nargs = 25;
  op_arg args[25];

  args[0] = arg0;
  arg1.idx = 0;
  args[1] = arg1;
  for ( int v=1; v<2; v++ ){
    args[1 + v] = op_arg_dat(arg1.dat, v, arg1.map, 60, "double", OP_READ);
  }

  arg3.idx = 0;
  args[3] = arg3;
  for ( int v=1; v<2; v++ ){
    args[3 + v] = op_arg_dat(arg3.dat, v, arg3.map, 60, "double", OP_READ);
  }

  arg5.idx = 0;
  args[5] = arg5;
  for ( int v=1; v<2; v++ ){
    args[5 + v] = op_arg_dat(arg5.dat, v, arg5.map, 60, "double", OP_READ);
  }

  arg7.idx = 0;
  args[7] = arg7;
  for ( int v=1; v<2; v++ ){
    args[7 + v] = op_arg_dat(arg7.dat, v, arg7.map, 60, "double", OP_READ);
  }

  arg9.idx = 0;
  args[9] = arg9;
  for ( int v=1; v<2; v++ ){
    args[9 + v] = op_arg_dat(arg9.dat, v, arg9.map, 60, "double", OP_READ);
  }

  arg11.idx = 0;
  args[11] = arg11;
  for ( int v=1; v<2; v++ ){
    args[11 + v] = op_arg_dat(arg11.dat, v, arg11.map, 60, "double", OP_READ);
  }

  arg13.idx = 0;
  args[13] = arg13;
  for ( int v=1; v<2; v++ ){
    args[13 + v] = op_arg_dat(arg13.dat, v, arg13.map, 60, "double", OP_INC);
  }

  arg15.idx = 0;
  args[15] = arg15;
  for ( int v=1; v<2; v++ ){
    args[15 + v] = op_arg_dat(arg15.dat, v, arg15.map, 60, "double", OP_INC);
  }

  arg17.idx = 0;
  args[17] = arg17;
  for ( int v=1; v<2; v++ ){
    args[17 + v] = op_arg_dat(arg17.dat, v, arg17.map, 60, "double", OP_INC);
  }

  arg19.idx = 0;
  args[19] = arg19;
  for ( int v=1; v<2; v++ ){
    args[19 + v] = op_arg_dat(arg19.dat, v, arg19.map, 60, "double", OP_INC);
  }

  arg21.idx = 0;
  args[21] = arg21;
  for ( int v=1; v<2; v++ ){
    args[21 + v] = op_arg_dat(arg21.dat, v, arg21.map, 60, "double", OP_INC);
  }

  arg23.idx = 0;
  args[23] = arg23;
  for ( int v=1; v<2; v++ ){
    args[23 + v] = op_arg_dat(arg23.dat, v, arg23.map, 60, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(11);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[11].name      = name;
  OP_kernels[11].count    += 1;

  int  ninds   = 12;
  int  inds[25] = {-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_grad_faces\n");
  }

  // get plan
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_11
    int part_size = OP_PART_SIZE_11;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_11
    int nthread = OP_BLOCK_SIZE_11;
  #else
    int nthread = OP_block_size;
  #endif


  int ncolors = 0;
  int set_size1 = set->size + set->exec_size;

  if (set_size >0) {

    //Set up typed device pointers for OpenMP
    int *map1 = arg1.map_data_d;
     int map1size = arg1.map->dim * set_size1;

    int* data0 = (int*)arg0.data_d;
    int dat0size = getSetSizeFromOpArg(&arg0) * arg0.dat->dim;
    double *data1 = (double *)arg1.data_d;
    int dat1size = getSetSizeFromOpArg(&arg1) * arg1.dat->dim;
    double *data3 = (double *)arg3.data_d;
    int dat3size = getSetSizeFromOpArg(&arg3) * arg3.dat->dim;
    double *data5 = (double *)arg5.data_d;
    int dat5size = getSetSizeFromOpArg(&arg5) * arg5.dat->dim;
    double *data7 = (double *)arg7.data_d;
    int dat7size = getSetSizeFromOpArg(&arg7) * arg7.dat->dim;
    double *data9 = (double *)arg9.data_d;
    int dat9size = getSetSizeFromOpArg(&arg9) * arg9.dat->dim;
    double *data11 = (double *)arg11.data_d;
    int dat11size = getSetSizeFromOpArg(&arg11) * arg11.dat->dim;
    double *data13 = (double *)arg13.data_d;
    int dat13size = getSetSizeFromOpArg(&arg13) * arg13.dat->dim;
    double *data15 = (double *)arg15.data_d;
    int dat15size = getSetSizeFromOpArg(&arg15) * arg15.dat->dim;
    double *data17 = (double *)arg17.data_d;
    int dat17size = getSetSizeFromOpArg(&arg17) * arg17.dat->dim;
    double *data19 = (double *)arg19.data_d;
    int dat19size = getSetSizeFromOpArg(&arg19) * arg19.dat->dim;
    double *data21 = (double *)arg21.data_d;
    int dat21size = getSetSizeFromOpArg(&arg21) * arg21.dat->dim;
    double *data23 = (double *)arg23.data_d;
    int dat23size = getSetSizeFromOpArg(&arg23) * arg23.dat->dim;

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

      gauss_grad_faces_omp4_kernel(
        data0,
        dat0size,
        map1,
        map1size,
        data1,
        dat1size,
        data3,
        dat3size,
        data5,
        dat5size,
        data7,
        dat7size,
        data9,
        dat9size,
        data11,
        dat11size,
        data13,
        dat13size,
        data15,
        dat15size,
        data17,
        dat17size,
        data19,
        dat19size,
        data21,
        dat21size,
        data23,
        dat23size,
        col_reord,
        set_size1,
        start,
        end,
        part_size!=0?(end-start-1)/part_size+1:(end-start-1)/nthread,
        nthread);

    }
    OP_kernels[11].transfer  += Plan->transfer;
    OP_kernels[11].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[11].time     += wall_t2 - wall_t1;
}
