//
// auto-generated by op2.py
//

//user function
//user function

void advection_bc_omp4_kernel(
  int *data0,
  int dat0size,
  int *data1,
  int dat1size,
  double *arg2,
  int *arg3,
  int *map4,
  int map4size,
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
  double *data9,
  int dat9size,
  double *data10,
  int dat10size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_advection_bc(char const *name, op_set set,
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
  op_arg arg10){

  double*arg2h = (double *)arg2.data;
  int*arg3h = (int *)arg3.data;
  int nargs = 11;
  op_arg args[11];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(34);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[34].name      = name;
  OP_kernels[34].count    += 1;

  int  ninds   = 7;
  int  inds[11] = {-1,-1,-1,-1,0,1,2,3,4,5,6};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: advection_bc\n");
  }

  // get plan
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_34
    int part_size = OP_PART_SIZE_34;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_34
    int nthread = OP_BLOCK_SIZE_34;
  #else
    int nthread = OP_block_size;
  #endif

  double arg2_l = arg2h[0];
  int arg3_l = arg3h[0];

  int ncolors = 0;
  int set_size1 = set->size + set->exec_size;

  if (set_size >0) {

    //Set up typed device pointers for OpenMP
    int *map4 = arg4.map_data_d;
     int map4size = arg4.map->dim * set_size1;

    int* data0 = (int*)arg0.data_d;
    int dat0size = getSetSizeFromOpArg(&arg0) * arg0.dat->dim;
    int* data1 = (int*)arg1.data_d;
    int dat1size = getSetSizeFromOpArg(&arg1) * arg1.dat->dim;
    double *data4 = (double *)arg4.data_d;
    int dat4size = getSetSizeFromOpArg(&arg4) * arg4.dat->dim;
    double *data5 = (double *)arg5.data_d;
    int dat5size = getSetSizeFromOpArg(&arg5) * arg5.dat->dim;
    double *data6 = (double *)arg6.data_d;
    int dat6size = getSetSizeFromOpArg(&arg6) * arg6.dat->dim;
    double *data7 = (double *)arg7.data_d;
    int dat7size = getSetSizeFromOpArg(&arg7) * arg7.dat->dim;
    double *data8 = (double *)arg8.data_d;
    int dat8size = getSetSizeFromOpArg(&arg8) * arg8.dat->dim;
    double *data9 = (double *)arg9.data_d;
    int dat9size = getSetSizeFromOpArg(&arg9) * arg9.dat->dim;
    double *data10 = (double *)arg10.data_d;
    int dat10size = getSetSizeFromOpArg(&arg10) * arg10.dat->dim;

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

      advection_bc_omp4_kernel(
        data0,
        dat0size,
        data1,
        dat1size,
        &arg2_l,
        &arg3_l,
        map4,
        map4size,
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
        data9,
        dat9size,
        data10,
        dat10size,
        col_reord,
        set_size1,
        start,
        end,
        part_size!=0?(end-start-1)/part_size+1:(end-start-1)/nthread,
        nthread);

    }
    OP_kernels[34].transfer  += Plan->transfer;
    OP_kernels[34].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[34].time     += wall_t2 - wall_t1;
}
