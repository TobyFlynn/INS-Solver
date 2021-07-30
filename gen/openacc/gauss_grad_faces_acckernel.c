//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void gauss_grad_faces_openacc( const int *edgeNum, const double **mDx0,
                             const double **mDy0, const double **mDx1,
                             const double **mDy1, const double **mDx2,
                             const double **mDy2, double **pDx0, double **pDy0,
                             double **pDx1, double **pDy1, double **pDx2,
                             double **pDy2) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  for(int m = 0; m < 4; m++) {
    for(int n = 0; n < 6; n++) {
      int indL = m * 6 + n;
      int indR = m * 6 + n;

      if(edgeL == 0) {
        if(edgeR == 0) {
          pDx0[0][indL] += mDx0[1][indR];
          pDy0[0][indL] += mDy0[1][indR];
          pDx0[1][indR] += mDx0[0][indL];
          pDy0[1][indR] += mDy0[0][indL];
        } else if(edgeR == 1) {
          pDx0[0][indL] += mDx1[1][indR];
          pDy0[0][indL] += mDy1[1][indR];
          pDx1[1][indR] += mDx0[0][indL];
          pDy1[1][indR] += mDy0[0][indL];
        } else {
          pDx0[0][indL] += mDx2[1][indR];
          pDy0[0][indL] += mDy2[1][indR];
          pDx2[1][indR] += mDx0[0][indL];
          pDy2[1][indR] += mDy0[0][indL];
        }
      } else if(edgeL == 1) {
        if(edgeR == 0) {
          pDx1[0][indL] += mDx0[1][indR];
          pDy1[0][indL] += mDy0[1][indR];
          pDx0[1][indR] += mDx1[0][indL];
          pDy0[1][indR] += mDy1[0][indL];
        } else if(edgeR == 1) {
          pDx1[0][indL] += mDx1[1][indR];
          pDy1[0][indL] += mDy1[1][indR];
          pDx1[1][indR] += mDx1[0][indL];
          pDy1[1][indR] += mDy1[0][indL];
        } else {
          pDx1[0][indL] += mDx2[1][indR];
          pDy1[0][indL] += mDy2[1][indR];
          pDx2[1][indR] += mDx1[0][indL];
          pDy2[1][indR] += mDy1[0][indL];
        }
      } else {
        if(edgeR == 0) {
          pDx2[0][indL] += mDx0[1][indR];
          pDy2[0][indL] += mDy0[1][indR];
          pDx0[1][indR] += mDx2[0][indL];
          pDy0[1][indR] += mDy2[0][indL];
        } else if(edgeR == 1) {
          pDx2[0][indL] += mDx1[1][indR];
          pDy2[0][indL] += mDy1[1][indR];
          pDx1[1][indR] += mDx2[0][indL];
          pDy1[1][indR] += mDy2[0][indL];
        } else {
          pDx2[0][indL] += mDx2[1][indR];
          pDy2[0][indL] += mDy2[1][indR];
          pDx2[1][indR] += mDx2[0][indL];
          pDy2[1][indR] += mDy2[0][indL];
        }
      }
    }
  }
}

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
    args[1 + v] = op_arg_dat(arg1.dat, v, arg1.map, 24, "double", OP_READ);
  }

  arg3.idx = 0;
  args[3] = arg3;
  for ( int v=1; v<2; v++ ){
    args[3 + v] = op_arg_dat(arg3.dat, v, arg3.map, 24, "double", OP_READ);
  }

  arg5.idx = 0;
  args[5] = arg5;
  for ( int v=1; v<2; v++ ){
    args[5 + v] = op_arg_dat(arg5.dat, v, arg5.map, 24, "double", OP_READ);
  }

  arg7.idx = 0;
  args[7] = arg7;
  for ( int v=1; v<2; v++ ){
    args[7 + v] = op_arg_dat(arg7.dat, v, arg7.map, 24, "double", OP_READ);
  }

  arg9.idx = 0;
  args[9] = arg9;
  for ( int v=1; v<2; v++ ){
    args[9 + v] = op_arg_dat(arg9.dat, v, arg9.map, 24, "double", OP_READ);
  }

  arg11.idx = 0;
  args[11] = arg11;
  for ( int v=1; v<2; v++ ){
    args[11 + v] = op_arg_dat(arg11.dat, v, arg11.map, 24, "double", OP_READ);
  }

  arg13.idx = 0;
  args[13] = arg13;
  for ( int v=1; v<2; v++ ){
    args[13 + v] = op_arg_dat(arg13.dat, v, arg13.map, 24, "double", OP_INC);
  }

  arg15.idx = 0;
  args[15] = arg15;
  for ( int v=1; v<2; v++ ){
    args[15 + v] = op_arg_dat(arg15.dat, v, arg15.map, 24, "double", OP_INC);
  }

  arg17.idx = 0;
  args[17] = arg17;
  for ( int v=1; v<2; v++ ){
    args[17 + v] = op_arg_dat(arg17.dat, v, arg17.map, 24, "double", OP_INC);
  }

  arg19.idx = 0;
  args[19] = arg19;
  for ( int v=1; v<2; v++ ){
    args[19 + v] = op_arg_dat(arg19.dat, v, arg19.map, 24, "double", OP_INC);
  }

  arg21.idx = 0;
  args[21] = arg21;
  for ( int v=1; v<2; v++ ){
    args[21 + v] = op_arg_dat(arg21.dat, v, arg21.map, 24, "double", OP_INC);
  }

  arg23.idx = 0;
  args[23] = arg23;
  for ( int v=1; v<2; v++ ){
    args[23 + v] = op_arg_dat(arg23.dat, v, arg23.map, 24, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(9);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[9].name      = name;
  OP_kernels[9].count    += 1;

  int  ninds   = 12;
  int  inds[25] = {-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_grad_faces\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_9
    int part_size = OP_PART_SIZE_9;
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
    double *data17 = (double *)arg17.data_d;
    double *data19 = (double *)arg19.data_d;
    double *data21 = (double *)arg21.data_d;
    double *data23 = (double *)arg23.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map1,data0,data1,data3,data5,data7,data9,data11,data13,data15,data17,data19,data21,data23)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map1idx;
        int map2idx;
        map1idx = map1[n + set_size1 * 0];
        map2idx = map1[n + set_size1 * 1];

        const double* arg1_vec[] = {
           &data1[24 * map1idx],
           &data1[24 * map2idx]};
        const double* arg3_vec[] = {
           &data3[24 * map1idx],
           &data3[24 * map2idx]};
        const double* arg5_vec[] = {
           &data5[24 * map1idx],
           &data5[24 * map2idx]};
        const double* arg7_vec[] = {
           &data7[24 * map1idx],
           &data7[24 * map2idx]};
        const double* arg9_vec[] = {
           &data9[24 * map1idx],
           &data9[24 * map2idx]};
        const double* arg11_vec[] = {
           &data11[24 * map1idx],
           &data11[24 * map2idx]};
        double* arg13_vec[] = {
           &data13[24 * map1idx],
           &data13[24 * map2idx]};
        double* arg15_vec[] = {
           &data15[24 * map1idx],
           &data15[24 * map2idx]};
        double* arg17_vec[] = {
           &data17[24 * map1idx],
           &data17[24 * map2idx]};
        double* arg19_vec[] = {
           &data19[24 * map1idx],
           &data19[24 * map2idx]};
        double* arg21_vec[] = {
           &data21[24 * map1idx],
           &data21[24 * map2idx]};
        double* arg23_vec[] = {
           &data23[24 * map1idx],
           &data23[24 * map2idx]};

        gauss_grad_faces_openacc(
          &data0[2 * n],
          arg1_vec,
          arg3_vec,
          arg5_vec,
          arg7_vec,
          arg9_vec,
          arg11_vec,
          arg13_vec,
          arg15_vec,
          arg17_vec,
          arg19_vec,
          arg21_vec,
          arg23_vec);
      }

    }
    OP_kernels[9].transfer  += Plan->transfer;
    OP_kernels[9].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[9].time     += wall_t2 - wall_t1;
}
