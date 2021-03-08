//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void set_tau_openacc( const int *edgeNum, const double **x, const double **y,
                    const double **J, const double **sJ, double **tau) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse;

  if(edgeR == 0) {
    if(edgeL == 0) {
      reverse = !(x[0][0] == x[1][0] && y[0][0] == y[1][0]);
    } else if(edgeL == 1) {
      reverse = !(x[0][1] == x[1][0] && y[0][1] == y[1][0]);
    } else {
      reverse = !(x[0][2] == x[1][0] && y[0][2] == y[1][0]);
    }
  } else if(edgeR == 1) {
    if(edgeL == 0) {
      reverse = !(x[0][0] == x[1][1] && y[0][0] == y[1][1]);
    } else if(edgeL == 1) {
      reverse = !(x[0][1] == x[1][1] && y[0][1] == y[1][1]);
    } else {
      reverse = !(x[0][2] == x[1][1] && y[0][2] == y[1][1]);
    }
  } else {
    if(edgeL == 0) {
      reverse = !(x[0][0] == x[1][2] && y[0][0] == y[1][2]);
    } else if(edgeL == 1) {
      reverse = !(x[0][1] == x[1][2] && y[0][1] == y[1][2]);
    } else {
      reverse = !(x[0][2] == x[1][2] && y[0][2] == y[1][2]);
    }
  }

  int exIndL = 0;
  if(edgeL == 1) exIndL = 5;
  else if(edgeL == 2) exIndL = 2 * 5;

  int exIndR = 0;
  if(edgeR == 1) exIndR = 5;
  else if(edgeR == 2) exIndR = 2 * 5;

  int *fmaskR;

  if(edgeR == 0) {
    fmaskR = FMASK;
  } else if(edgeR == 1) {
    fmaskR = &FMASK[5];
  } else {
    fmaskR = &FMASK[2 * 5];
  }

  int *fmaskL;

  if(edgeL == 0) {
    fmaskL = FMASK;
  } else if(edgeL == 1) {
    fmaskL = &FMASK[5];
  } else {
    fmaskL = &FMASK[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int rIndF, lIndF, rInd, lInd;
    if(reverse) {
      rIndF = fmaskR[5 - i - 1];
      rInd = exIndR + 5 - i - 1;
    } else {
      rIndF = fmaskR[i];
      rInd = exIndR + i;
    }
    lIndF = fmaskL[i];
    lInd = exIndL + i;

    double lH = 2.0 * J[0][lIndF] / sJ[0][lInd];
    double rH = 2.0 * J[1][rIndF] / sJ[1][rInd];
    if(lH < rH) {
      tau[0][lInd] += 15.0 / lH;
      tau[1][rInd] += 15.0 / lH;
    } else {
      tau[0][lInd] += 15.0 / rH;
      tau[1][rInd] += 15.0 / rH;
    }
  }
}

// host stub function
void op_par_loop_set_tau(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3,
  op_arg arg5,
  op_arg arg7,
  op_arg arg9){

  int nargs = 11;
  op_arg args[11];

  args[0] = arg0;
  arg1.idx = 0;
  args[1] = arg1;
  for ( int v=1; v<2; v++ ){
    args[1 + v] = op_arg_dat(arg1.dat, v, arg1.map, 3, "double", OP_READ);
  }

  arg3.idx = 0;
  args[3] = arg3;
  for ( int v=1; v<2; v++ ){
    args[3 + v] = op_arg_dat(arg3.dat, v, arg3.map, 3, "double", OP_READ);
  }

  arg5.idx = 0;
  args[5] = arg5;
  for ( int v=1; v<2; v++ ){
    args[5 + v] = op_arg_dat(arg5.dat, v, arg5.map, 15, "double", OP_READ);
  }

  arg7.idx = 0;
  args[7] = arg7;
  for ( int v=1; v<2; v++ ){
    args[7 + v] = op_arg_dat(arg7.dat, v, arg7.map, 15, "double", OP_READ);
  }

  arg9.idx = 0;
  args[9] = arg9;
  for ( int v=1; v<2; v++ ){
    args[9 + v] = op_arg_dat(arg9.dat, v, arg9.map, 15, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(25);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[25].name      = name;
  OP_kernels[25].count    += 1;

  int  ninds   = 5;
  int  inds[11] = {-1,0,0,1,1,2,2,3,3,4,4};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: set_tau\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_25
    int part_size = OP_PART_SIZE_25;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map1 = arg1.map_data_d;

    int* data0 = (int*)arg0.data_d;
    double *data1 = (double *)arg1.data_d;
    double *data3 = (double *)arg3.data_d;
    double *data5 = (double *)arg5.data_d;
    double *data7 = (double *)arg7.data_d;
    double *data9 = (double *)arg9.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map1,data0,data1,data3,data5,data7,data9)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map1idx;
        int map2idx;
        map1idx = map1[n + set_size1 * 0];
        map2idx = map1[n + set_size1 * 1];

        const double* arg1_vec[] = {
           &data1[3 * map1idx],
           &data1[3 * map2idx]};
        const double* arg3_vec[] = {
           &data3[3 * map1idx],
           &data3[3 * map2idx]};
        const double* arg5_vec[] = {
           &data5[15 * map1idx],
           &data5[15 * map2idx]};
        const double* arg7_vec[] = {
           &data7[15 * map1idx],
           &data7[15 * map2idx]};
        double* arg9_vec[] = {
           &data9[15 * map1idx],
           &data9[15 * map2idx]};

        set_tau_openacc(
          &data0[2 * n],
          arg1_vec,
          arg3_vec,
          arg5_vec,
          arg7_vec,
          arg9_vec);
      }

    }
    OP_kernels[25].transfer  += Plan->transfer;
    OP_kernels[25].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[25].time     += wall_t2 - wall_t1;
}
