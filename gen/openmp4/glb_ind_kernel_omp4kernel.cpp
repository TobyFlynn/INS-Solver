//
// auto-generated by op2.py
//

//user function
//user function

void glb_ind_kernel_omp4_kernel(
  int *map0,
  int map0size,
  int *data2,
  int dat2size,
  int *data3,
  int dat3size,
  int *data0,
  int dat0size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_glb_ind_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  arg0.idx = 0;
  args[0] = arg0;
  for ( int v=1; v<2; v++ ){
    args[0 + v] = op_arg_dat(arg0.dat, v, arg0.map, 1, "int", OP_READ);
  }

  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(15);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[15].name      = name;
  OP_kernels[15].count    += 1;

  int  ninds   = 1;
  int  inds[4] = {0,0,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: glb_ind_kernel\n");
  }

  // get plan
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_15
    int part_size = OP_PART_SIZE_15;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_15
    int nthread = OP_BLOCK_SIZE_15;
  #else
    int nthread = OP_block_size;
  #endif


  int ncolors = 0;
  int set_size1 = set->size + set->exec_size;

  if (set_size >0) {

    //Set up typed device pointers for OpenMP
    int *map0 = arg0.map_data_d;
     int map0size = arg0.map->dim * set_size1;

    int* data2 = (int*)arg2.data_d;
    int dat2size = getSetSizeFromOpArg(&arg2) * arg2.dat->dim;
    int* data3 = (int*)arg3.data_d;
    int dat3size = getSetSizeFromOpArg(&arg3) * arg3.dat->dim;
    int *data0 = (int *)arg0.data_d;
    int dat0size = getSetSizeFromOpArg(&arg0) * arg0.dat->dim;

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

      glb_ind_kernel_omp4_kernel(
        map0,
        map0size,
        data2,
        dat2size,
        data3,
        dat3size,
        data0,
        dat0size,
        col_reord,
        set_size1,
        start,
        end,
        part_size!=0?(end-start-1)/part_size+1:(end-start-1)/nthread,
        nthread);

    }
    OP_kernels[15].transfer  += Plan->transfer;
    OP_kernels[15].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[15].time     += wall_t2 - wall_t1;
}
