//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void init_gauss_grad3_2_openacc( const int *edgeNum,
                             const double *nxL, const double *nxR,
                             const double *nyL, const double *nyR,
                             const double *Dx0L, const double *Dx0R,
                             const double *Dy0L, const double *Dy0R,
                             const double *Dx1L, const double *Dx1R,
                             const double *Dy1L, const double *Dy1R,
                             const double *Dx2L, const double *Dx2R,
                             const double *Dy2L, const double *Dy2R,
                             const double *factL, const double *factR,
                             double *dL, double *dR) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  const double *DxL, *DyL, *DxR, *DyR;

  if(edgeL == 0) {
    DxL = Dx0L;
    DyL = Dy0L;
  } else if(edgeL == 1) {
    DxL = Dx1L;
    DyL = Dy1L;
  } else {
    DxL = Dx2L;
    DyL = Dy2L;
  }

  if(edgeR == 0) {
    DxR = Dx0R;
    DyR = Dy0R;
  } else if(edgeR == 1) {
    DxR = Dx1R;
    DyR = Dy1R;
  } else {
    DxR = Dx2R;
    DyR = Dy2R;
  }

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int ind  = m * 10 + n;
      int indL = edgeL * 6 + m;
      int indR = edgeR * 6 + m;

      dL[ind] = nxL[indL] * factL[indL] * DxL[ind] + nyL[indL] * factL[indL] * DyL[ind];
      dR[ind] = nxR[indR] * factR[indR] * DxR[ind] + nyR[indR] * factR[indR] * DyR[ind];
    }
  }
}

// host stub function
void op_par_loop_init_gauss_grad3_2(char const *name, op_set set,
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
  op_arg arg18,
  op_arg arg19,
  op_arg arg20){

  int nargs = 21;
  op_arg args[21];

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
  args[19] = arg19;
  args[20] = arg20;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(16);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[16].name      = name;
  OP_kernels[16].count    += 1;

  int  ninds   = 9;
  int  inds[21] = {-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: init_gauss_grad3_2\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_16
    int part_size = OP_PART_SIZE_16;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map1 = arg1.map_data_d;

    int* data0 = (int*)arg0.data_d;
    double* data19 = (double*)arg19.data_d;
    double* data20 = (double*)arg20.data_d;
    double *data1 = (double *)arg1.data_d;
    double *data3 = (double *)arg3.data_d;
    double *data5 = (double *)arg5.data_d;
    double *data7 = (double *)arg7.data_d;
    double *data9 = (double *)arg9.data_d;
    double *data11 = (double *)arg11.data_d;
    double *data13 = (double *)arg13.data_d;
    double *data15 = (double *)arg15.data_d;
    double *data17 = (double *)arg17.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map1,data0,data19,data20,data1,data3,data5,data7,data9,data11,data13,data15,data17)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map1idx;
        int map2idx;
        map1idx = map1[n + set_size1 * 0];
        map2idx = map1[n + set_size1 * 1];


        init_gauss_grad3_2_openacc(
          &data0[2 * n],
          &data1[18 * map1idx],
          &data1[18 * map2idx],
          &data3[18 * map1idx],
          &data3[18 * map2idx],
          &data5[60 * map1idx],
          &data5[60 * map2idx],
          &data7[60 * map1idx],
          &data7[60 * map2idx],
          &data9[60 * map1idx],
          &data9[60 * map2idx],
          &data11[60 * map1idx],
          &data11[60 * map2idx],
          &data13[60 * map1idx],
          &data13[60 * map2idx],
          &data15[60 * map1idx],
          &data15[60 * map2idx],
          &data17[18 * map1idx],
          &data17[18 * map2idx],
          &data19[60 * n],
          &data20[60 * n]);
      }

    }
    OP_kernels[16].transfer  += Plan->transfer;
    OP_kernels[16].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[16].time     += wall_t2 - wall_t1;
}
