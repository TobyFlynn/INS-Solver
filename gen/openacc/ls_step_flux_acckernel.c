//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void ls_step_flux_openacc( const int *edgeNum, const bool *rev,
                         const double **fscale, const double **step,
                         double **flux) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exInd = edgeL * 4;
  int *fmaskL = &FMASK[edgeL * 4];
  int *fmaskR = &FMASK[edgeR * 4];

  for(int i = 0; i < 4; i++) {
    int lInd = fmaskL[i];
    int rInd;
    if(reverse) {
      rInd = fmaskR[4 - i - 1];
    } else {
      rInd = fmaskR[i];
    }
    double tmp = step[0][lInd] - (step[0][lInd] + step[1][rInd]) / 2.0;
    flux[0][exInd + i] += fscale[0][exInd + i] * tmp;
  }

  exInd = edgeR * 4;

  for(int i = 0; i < 4; i++) {
    int rInd = fmaskR[i];
    int lInd;
    if(reverse) {
      lInd = fmaskL[4 - i - 1];
    } else {
      lInd = fmaskL[i];
    }
    double tmp = step[1][rInd] - (step[0][lInd] + step[1][rInd]) / 2.0;
    flux[1][exInd + i] += fscale[1][exInd + i] * tmp;
  }
}

// host stub function
void op_par_loop_ls_step_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6){

  int nargs = 8;
  op_arg args[8];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 12, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 10, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 12, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(35);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[35].name      = name;
  OP_kernels[35].count    += 1;

  int  ninds   = 3;
  int  inds[8] = {-1,-1,0,0,1,1,2,2};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_step_flux\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_35
    int part_size = OP_PART_SIZE_35;
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
    double *data6 = (double *)arg6.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data4,data6)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        int map3idx;
        map2idx = map2[n + set_size1 * 0];
        map3idx = map2[n + set_size1 * 1];

        const double* arg2_vec[] = {
           &data2[12 * map2idx],
           &data2[12 * map3idx]};
        const double* arg4_vec[] = {
           &data4[10 * map2idx],
           &data4[10 * map3idx]};
        double* arg6_vec[] = {
           &data6[12 * map2idx],
           &data6[12 * map3idx]};

        ls_step_flux_openacc(
          &data0[2 * n],
          &data1[1 * n],
          arg2_vec,
          arg4_vec,
          arg6_vec);
      }

    }
    OP_kernels[35].transfer  += Plan->transfer;
    OP_kernels[35].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[35].time     += wall_t2 - wall_t1;
}
