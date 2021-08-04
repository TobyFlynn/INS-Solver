//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void ls_advec_edges_openacc( const int *edgeNum, const bool *rev,
                           const double **q, double **exQ) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exInd = edgeL * 4;
  int *fmask = &FMASK[edgeR * 4];

  for(int i = 0; i < 4; i++) {
    int rInd;
    if(reverse) {
      rInd = fmask[4 - i - 1];
    } else {
      rInd = fmask[i];
    }
    exQ[0][exInd + i] += q[1][rInd];
  }

  exInd = edgeR * 4;
  fmask = &FMASK[edgeL * 4];

  for(int i = 0; i < 4; i++) {
    int lInd;
    if(reverse) {
      lInd = fmask[4 - i - 1];
    } else {
      lInd = fmask[i];
    }
    exQ[1][exInd + i] += q[0][lInd];
  }
}

// host stub function
void op_par_loop_ls_advec_edges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 10, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 12, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(51);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[51].name      = name;
  OP_kernels[51].count    += 1;

  int  ninds   = 2;
  int  inds[6] = {-1,-1,0,0,1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_advec_edges\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_51
    int part_size = OP_PART_SIZE_51;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map2 = arg2.map_data_d;

    int* data0 = (int*)arg0.data_d;
    bool* data1 = (bool*)arg1.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data4 = (double *)arg4.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data4)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        int map3idx;
        map2idx = map2[n + set_size1 * 0];
        map3idx = map2[n + set_size1 * 1];

        const double* arg2_vec[] = {
           &data2[10 * map2idx],
           &data2[10 * map3idx]};
        double* arg4_vec[] = {
           &data4[12 * map2idx],
           &data4[12 * map3idx]};

        ls_advec_edges_openacc(
          &data0[2 * n],
          &data1[1 * n],
          arg2_vec,
          arg4_vec);
      }

    }
    OP_kernels[51].transfer  += Plan->transfer;
    OP_kernels[51].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[51].time     += wall_t2 - wall_t1;
}
