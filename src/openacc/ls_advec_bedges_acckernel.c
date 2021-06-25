//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void ls_advec_bedges_openacc( const int *bedge_type, const int *bedgeNum,
                            const double *x, const double *y, const double *q,
                            double *exQ) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 5;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 5;
  }

  int *fmask;

  if(*bedgeNum == 0) {
    fmask = FMASK;
  } else if(*bedgeNum == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    exQ[exInd + i] += q[fmask[i]];
  }
}

// host stub function
void op_par_loop_ls_advec_bedges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

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
  op_timing_realloc(52);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[52].name      = name;
  OP_kernels[52].count    += 1;

  int  ninds   = 4;
  int  inds[6] = {-1,-1,0,1,2,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_advec_bedges\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_52
    int part_size = OP_PART_SIZE_52;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map2 = arg2.map_data_d;

    int* data0 = (int*)arg0.data_d;
    int* data1 = (int*)arg1.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data3 = (double *)arg3.data_d;
    double *data4 = (double *)arg4.data_d;
    double *data5 = (double *)arg5.data_d;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;
    int set_size1 = set->size + set->exec_size;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data3,data4,data5)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        map2idx = map2[n + set_size1 * 0];


        ls_advec_bedges_openacc(
          &data0[1 * n],
          &data1[1 * n],
          &data2[15 * map2idx],
          &data3[15 * map2idx],
          &data4[15 * map2idx],
          &data5[15 * map2idx]);
      }

    }
    OP_kernels[52].transfer  += Plan->transfer;
    OP_kernels[52].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[52].time     += wall_t2 - wall_t1;
}