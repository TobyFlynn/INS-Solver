//
// auto-generated by op2.py
//

//user function
//user function

void advection_faces_omp4_kernel(
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
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_advection_faces(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6,
  op_arg arg8){

  int nargs = 10;
  op_arg args[10];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 15, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 15, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 15, "double", OP_INC);
  }

  arg8.idx = 0;
  args[8] = arg8;
  for ( int v=1; v<2; v++ ){
    args[8 + v] = op_arg_dat(arg8.dat, v, arg8.map, 15, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(39);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[39].name      = name;
  OP_kernels[39].count    += 1;

  int  ninds   = 4;
  int  inds[10] = {-1,-1,0,0,1,1,2,2,3,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: advection_faces\n");
  }

  // get plan
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_39
    int part_size = OP_PART_SIZE_39;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_39
    int nthread = OP_BLOCK_SIZE_39;
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

      advection_faces_omp4_kernel(
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
        col_reord,
        set_size1,
        start,
        end,
        part_size!=0?(end-start-1)/part_size+1:(end-start-1)/nthread,
        nthread);

    }
    OP_kernels[39].transfer  += Plan->transfer;
    OP_kernels[39].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[39].time     += wall_t2 - wall_t1;
}
