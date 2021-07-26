//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_op2_openacc( const int *edgeNum, const bool *rev,
                        const double *mD0L, const double *mD0R,
                        const double *mD1L, const double *mD1R,
                        const double *mD2L, const double *mD2R,
                        const double *pD0L, const double *pD0R,
                        const double *pD1L, const double *pD1R,
                        const double *pD2L, const double *pD2R,
                        const double *gVP0L, const double *gVP0R,
                        const double *gVP1L, const double *gVP1R,
                        const double *gVP2L, const double *gVP2R,
                        const double *sJL, const double *sJR,
                        const double *hL, const double *hR,
                        const double *gFactorL, const double *gFactorR,
                        const double *factorL, const double *factorR,
                        double *op1L, double *op1R,
                        double *op2L, double *op2R) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;


  const double *mDL, *mDR, *pDL, *pDR, *gVML, *gVMR, *gVPL, *gVPR;
  if(edgeL == 0) {
    mDL  = mD0L;
    pDL  = pD0L;
    gVML = gFInterp0_g;
    gVPL = gVP0L;
  } else if(edgeL == 1) {
    mDL  = mD1L;
    pDL  = pD1L;
    gVML = gFInterp1_g;
    gVPL = gVP1L;
  } else {
    mDL  = mD2L;
    pDL  = pD2L;
    gVML = gFInterp2_g;
    gVPL = gVP2L;
  }

  if(edgeR == 0) {
    mDR  = mD0R;
    pDR  = pD0R;
    gVMR = gFInterp0_g;
    gVPR = gVP0R;
  } else if(edgeR == 1) {
    mDR  = mD1R;
    pDR  = pD1R;
    gVMR = gFInterp1_g;
    gVPR = gVP1R;
  } else {
    mDR  = mD2R;
    pDR  = pD2R;
    gVMR = gFInterp2_g;
    gVPR = gVP2R;
  }



  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      op2L[c_ind] = 0.0;
      op2R[c_ind] = 0.0;
      for(int k = 0; k < 7; k++) {

        int b_ind = k * 15 + j;

        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);


        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;
        int factors_indLR;
        int factors_indRR;
        if(reverse) {
          factors_indLR = edgeL * 7 + 6 - k;
          factors_indRR = edgeR * 7 + 6 - k;
        } else {
          factors_indLR = edgeL * 7 + k;
          factors_indRR = edgeR * 7 + k;
        }

        op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * gFactorL[factors_indL] * mDL[b_ind];
        op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * gFactorR[factors_indR] * mDR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * gFactorR[factors_indRR] * pDL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * gFactorL[factors_indLR] * pDR[b_ind];
      }
    }
  }



  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      for(int k = 0; k < 7; k++) {

        int b_ind = k * 15 + j;

        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);

        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;

        op1L[c_ind] += -factorL[i] * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVML[b_ind];
        op1R[c_ind] += -factorR[i] * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVMR[b_ind];

        op2L[c_ind] += factorL[i] * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVPL[b_ind];
        op2R[c_ind] += factorR[i] * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVPR[b_ind];









      }
    }
  }

  double tauL[7];
  double tauR[7];
  double maxL = 0.0;
  double maxR = 0.0;
  for(int i = 0; i < 7; i++) {
    int indL = edgeL * 7 + i;
    int indR;
    if(reverse)
      indR = edgeR * 7 + 6 - i;
    else
      indR = edgeR * 7 + i;
    tauL[i] = 100 * 0.5 * 5 * 6 * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);

    if(maxL < tauL[i]) {
      maxL = tauL[i];
    }
  }
  for(int i = 0; i < 7; i++) {
    int indL;
    int indR = edgeR * 7 + i;
    if(reverse)
      indL = edgeL * 7 + 6 - i;
    else
      indL = edgeL * 7 + i;
    tauR[i] = 100 * 0.5 * 5 * 6 * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);

    if(maxR < tauR[i]) {
      maxR = tauR[i];
    }
  }

  for(int i = 0; i < 7; i++) {
    tauL[i] = maxL;
    tauR[i] = maxR;
  }



  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      for(int k = 0; k < 7; k++) {

        int b_ind = k * 15 + j;

        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);

        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;

        op1L[c_ind] += gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tauL[k] * gVML[b_ind];
        op1R[c_ind] += gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * tauR[k] * gVMR[b_ind];

        op2L[c_ind] += -gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tauL[k] * gVPL[b_ind];
        op2R[c_ind] += -gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * tauR[k] * gVPR[b_ind];
      }
    }
  }
}

// host stub function
void op_par_loop_poisson_op2(char const *name, op_set set,
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
  op_arg arg20,
  op_arg arg21,
  op_arg arg22,
  op_arg arg23,
  op_arg arg24,
  op_arg arg25,
  op_arg arg26,
  op_arg arg27,
  op_arg arg28,
  op_arg arg29,
  op_arg arg30,
  op_arg arg31){

  int nargs = 32;
  op_arg args[32];

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
  args[21] = arg21;
  args[22] = arg22;
  args[23] = arg23;
  args[24] = arg24;
  args[25] = arg25;
  args[26] = arg26;
  args[27] = arg27;
  args[28] = arg28;
  args[29] = arg29;
  args[30] = arg30;
  args[31] = arg31;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(20);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[20].name      = name;
  OP_kernels[20].count    += 1;

  int  ninds   = 14;
  int  inds[32] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_op2\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_20
    int part_size = OP_PART_SIZE_20;
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
    double* data30 = (double*)arg30.data_d;
    double* data31 = (double*)arg31.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data4 = (double *)arg4.data_d;
    double *data6 = (double *)arg6.data_d;
    double *data8 = (double *)arg8.data_d;
    double *data10 = (double *)arg10.data_d;
    double *data12 = (double *)arg12.data_d;
    double *data14 = (double *)arg14.data_d;
    double *data16 = (double *)arg16.data_d;
    double *data18 = (double *)arg18.data_d;
    double *data20 = (double *)arg20.data_d;
    double *data22 = (double *)arg22.data_d;
    double *data24 = (double *)arg24.data_d;
    double *data26 = (double *)arg26.data_d;
    double *data28 = (double *)arg28.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data30,data31,data2,data4,data6,data8,data10,data12,data14,data16,data18,data20,data22,data24,data26,data28)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        int map3idx;
        map2idx = map2[n + set_size1 * 0];
        map3idx = map2[n + set_size1 * 1];


        poisson_op2_openacc(
          &data0[2 * n],
          &data1[1 * n],
          &data2[105 * map2idx],
          &data2[105 * map3idx],
          &data4[105 * map2idx],
          &data4[105 * map3idx],
          &data6[105 * map2idx],
          &data6[105 * map3idx],
          &data8[105 * map2idx],
          &data8[105 * map3idx],
          &data10[105 * map2idx],
          &data10[105 * map3idx],
          &data12[105 * map2idx],
          &data12[105 * map3idx],
          &data14[105 * map2idx],
          &data14[105 * map3idx],
          &data16[105 * map2idx],
          &data16[105 * map3idx],
          &data18[105 * map2idx],
          &data18[105 * map3idx],
          &data20[21 * map2idx],
          &data20[21 * map3idx],
          &data22[1 * map2idx],
          &data22[1 * map3idx],
          &data24[21 * map2idx],
          &data24[21 * map3idx],
          &data26[15 * map2idx],
          &data26[15 * map3idx],
          &data28[225 * map2idx],
          &data28[225 * map3idx],
          &data30[225 * n],
          &data31[225 * n]);
      }

    }
    OP_kernels[20].transfer  += Plan->transfer;
    OP_kernels[20].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[20].time     += wall_t2 - wall_t1;
}
