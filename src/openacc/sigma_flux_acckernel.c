//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void sigma_flux_openacc( const int *edgeNum, const double **x, const double **y,
                       const double **sJ, const double **nx, const double **ny,
                       const double **s, double **sigFx, double **sigFy) {

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
  if(edgeL == 1) exIndL = 7;
  else if(edgeL == 2) exIndL = 2 * 7;

  int exIndR = 0;
  if(edgeR == 1) exIndR = 7;
  else if(edgeR == 2) exIndR = 2 * 7;

  for(int i = 0; i < 7; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + 7 - i - 1;
    } else {
      rInd = exIndR + i;
    }
    double flux = (s[0][lInd] + s[1][rInd]) / 2.0;
    sigFx[0][lInd] += gaussW_g[i] * sJ[0][lInd] * nx[0][lInd] * flux;
    sigFy[0][lInd] += gaussW_g[i] * sJ[0][lInd] * ny[0][lInd] * flux;
  }

  for(int i = 0; i < 7; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 7 - i - 1;
    } else {
      lInd = exIndL + i;
    }
    double flux = (s[0][lInd] + s[1][rInd]) / 2.0;
    sigFx[1][rInd] += gaussW_g[i] * sJ[1][rInd] * nx[1][rInd] * flux;
    sigFy[1][rInd] += gaussW_g[i] * sJ[1][rInd] * ny[1][rInd] * flux;
  }
}

// host stub function
void op_par_loop_sigma_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3,
  op_arg arg5,
  op_arg arg7,
  op_arg arg9,
  op_arg arg11,
  op_arg arg13,
  op_arg arg15){

  int nargs = 17;
  op_arg args[17];

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
    args[5 + v] = op_arg_dat(arg5.dat, v, arg5.map, 21, "double", OP_READ);
  }

  arg7.idx = 0;
  args[7] = arg7;
  for ( int v=1; v<2; v++ ){
    args[7 + v] = op_arg_dat(arg7.dat, v, arg7.map, 21, "double", OP_READ);
  }

  arg9.idx = 0;
  args[9] = arg9;
  for ( int v=1; v<2; v++ ){
    args[9 + v] = op_arg_dat(arg9.dat, v, arg9.map, 21, "double", OP_READ);
  }

  arg11.idx = 0;
  args[11] = arg11;
  for ( int v=1; v<2; v++ ){
    args[11 + v] = op_arg_dat(arg11.dat, v, arg11.map, 21, "double", OP_READ);
  }

  arg13.idx = 0;
  args[13] = arg13;
  for ( int v=1; v<2; v++ ){
    args[13 + v] = op_arg_dat(arg13.dat, v, arg13.map, 21, "double", OP_INC);
  }

  arg15.idx = 0;
  args[15] = arg15;
  for ( int v=1; v<2; v++ ){
    args[15 + v] = op_arg_dat(arg15.dat, v, arg15.map, 21, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(65);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[65].name      = name;
  OP_kernels[65].count    += 1;

  int  ninds   = 8;
  int  inds[17] = {-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: sigma_flux\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_65
    int part_size = OP_PART_SIZE_65;
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
    double *data11 = (double *)arg11.data_d;
    double *data13 = (double *)arg13.data_d;
    double *data15 = (double *)arg15.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map1,data0,data1,data3,data5,data7,data9,data11,data13,data15)
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
           &data5[21 * map1idx],
           &data5[21 * map2idx]};
        const double* arg7_vec[] = {
           &data7[21 * map1idx],
           &data7[21 * map2idx]};
        const double* arg9_vec[] = {
           &data9[21 * map1idx],
           &data9[21 * map2idx]};
        const double* arg11_vec[] = {
           &data11[21 * map1idx],
           &data11[21 * map2idx]};
        double* arg13_vec[] = {
           &data13[21 * map1idx],
           &data13[21 * map2idx]};
        double* arg15_vec[] = {
           &data15[21 * map1idx],
           &data15[21 * map2idx]};

        sigma_flux_openacc(
          &data0[2 * n],
          arg1_vec,
          arg3_vec,
          arg5_vec,
          arg7_vec,
          arg9_vec,
          arg11_vec,
          arg13_vec,
          arg15_vec);
      }

    }
    OP_kernels[65].transfer  += Plan->transfer;
    OP_kernels[65].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[65].time     += wall_t2 - wall_t1;
}
