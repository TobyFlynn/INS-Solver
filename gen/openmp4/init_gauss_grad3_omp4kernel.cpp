//
// auto-generated by op2.py
//

//user function
//user function

void init_gauss_grad3_omp4_kernel(
  int *data0,
  int dat0size,
  int *map1,
  int map1size,
  double *data17,
  int dat17size,
  double *data18,
  int dat18size,
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
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread);

// host stub function
void op_par_loop_init_gauss_grad3(char const *name, op_set set,
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
  op_arg arg15,
  op_arg arg16,
  op_arg arg17,
  op_arg arg18){

  int nargs = 19;
  op_arg args[19];

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
  args[16] = arg16;
  args[17] = arg17;
  args[18] = arg18;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(4);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[4].name      = name;
  OP_kernels[4].count    += 1;

  int  ninds   = 8;
  int  inds[19] = {-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: init_gauss_grad3\n");
  }

  // get plan
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_4
    int part_size = OP_PART_SIZE_4;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_4
    int nthread = OP_BLOCK_SIZE_4;
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
    double* data17 = (double*)arg17.data_d;
    int dat17size = getSetSizeFromOpArg(&arg17) * arg17.dat->dim;
    double* data18 = (double*)arg18.data_d;
    int dat18size = getSetSizeFromOpArg(&arg18) * arg18.dat->dim;
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

      init_gauss_grad3_omp4_kernel(
        data0,
        dat0size,
        map1,
        map1size,
        data17,
        dat17size,
        data18,
        dat18size,
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
        col_reord,
        set_size1,
        start,
        end,
        part_size!=0?(end-start-1)/part_size+1:(end-start-1)/nthread,
        nthread);

    }
    OP_kernels[4].transfer  += Plan->transfer;
    OP_kernels[4].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[4].time     += wall_t2 - wall_t1;
}
