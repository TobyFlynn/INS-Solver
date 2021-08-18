//
// auto-generated by op2.py
//

//user function
//user function
//#pragma acc routine
inline void poisson_op2_openacc( const int *edgeNum, const bool *rev,
                        const double *mDL, const double *mDR,
                        const double *pDL, const double *pDR,
                        const double *gVPL, const double *gVPR,
                        const double *sJL, const double *sJR,
                        const double *hL, const double *hR,
                        const double *gFactorL, const double *gFactorR,
                        const double *factorL, const double *factorR,
                        double *op1L, double *op1R,
                        double *op2L, double *op2R) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  const double *gVML, *gVMR;
  if(edgeL == 0) {
    gVML = gFInterp0_g;
  } else if(edgeL == 1) {
    gVML = gFInterp1_g;
  } else {
    gVML = gFInterp2_g;
  }

  if(edgeR == 0) {
    gVMR = gFInterp0_g;
  } else if(edgeR == 1) {
    gVMR = gFInterp1_g;
  } else {
    gVMR = gFInterp2_g;
  }



  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      op2L[c_ind] = 0.0;
      op2R[c_ind] = 0.0;
      for(int k = 0; k < 6; k++) {

        int b_ind = k * 10 + j;

        int ind = i * 6 + k;
        int a_ind = ((ind * 10) % (10 * 6)) + (ind / 6);


        int factors_indL = edgeL * 6 + k;
        int factors_indR = edgeR * 6 + k;
        int factors_indLR;
        int factors_indRR;
        if(reverse) {
          factors_indLR = edgeL * 6 + 6 - 1 - k;
          factors_indRR = edgeR * 6 + 6 - 1 - k;
        } else {
          factors_indLR = edgeL * 6 + k;
          factors_indRR = edgeR * 6 + k;
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



  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      for(int k = 0; k < 6; k++) {

        int b_ind = k * 10 + j;

        int ind = i * 6 + k;
        int a_ind = ((ind * 10) % (10 * 6)) + (ind / 6);

        int factors_indL = edgeL * 6 + k;
        int factors_indR = edgeR * 6 + k;
        int factors_indLR;
        int factors_indRR;
        if(reverse) {
          factors_indLR = edgeL * 6 + 6 - 1 - k;
          factors_indRR = edgeR * 6 + 6 - 1 - k;
        } else {
          factors_indLR = edgeL * 6 + k;
          factors_indRR = edgeR * 6 + k;
        }



















        op1L[c_ind] += -0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVML[b_ind];
        op1R[c_ind] += -0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVMR[b_ind];

        op2L[c_ind] += 0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g[k]
                       * sJL[factors_indL] * gVPL[b_ind];
        op2R[c_ind] += 0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g[k]
                       * sJR[factors_indR] * gVPR[b_ind];
      }
    }
  }

  double tauL[6];
  double tauR[6];
  double maxL = 0.0;
  double maxR = 0.0;
  for(int i = 0; i < 6; i++) {
    int indL = edgeL * 6 + i;
    int indR;
    if(reverse)
      indR = edgeR * 6 + 6 - 1 - i;
    else
      indR = edgeR * 6 + i;

    tauL[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);



  }
  for(int i = 0; i < 6; i++) {
    int indL;
    int indR = edgeR * 6 + i;
    if(reverse)
      indL = edgeL * 6 + 6 - 1 - i;
    else
      indL = edgeL * 6 + i;

    tauR[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(*hL * gFactorL[indL], *hR * gFactorR[indR]);



  }







  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      for(int k = 0; k < 6; k++) {

        int b_ind = k * 10 + j;

        int ind = i * 6 + k;
        int a_ind = ((ind * 10) % (10 * 6)) + (ind / 6);

        int factors_indL = edgeL * 6 + k;
        int factors_indR = edgeR * 6 + k;










        op1L[c_ind] += 0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tauL[k] * gVML[b_ind];
        op1R[c_ind] += 0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
                       * tauR[k] * gVMR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g[k] * sJL[factors_indL]
                       * tauL[k] * gVPL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g[k] * sJR[factors_indR]
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
  op_arg arg19){

  int nargs = 20;
  op_arg args[20];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(23);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[23].name      = name;
  OP_kernels[23].count    += 1;

  int  ninds   = 5;
  int  inds[20] = {-1,-1,-1,-1,-1,-1,-1,-1,0,0,1,1,2,2,3,3,4,4,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_op2\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_23
    int part_size = OP_PART_SIZE_23;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {


    //Set up typed device pointers for OpenACC
    int *map8 = arg8.map_data_d;

    int* data0 = (int*)arg0.data_d;
    bool* data1 = (bool*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    double* data4 = (double*)arg4.data_d;
    double* data5 = (double*)arg5.data_d;
    double* data6 = (double*)arg6.data_d;
    double* data7 = (double*)arg7.data_d;
    double* data18 = (double*)arg18.data_d;
    double* data19 = (double*)arg19.data_d;
    double *data8 = (double *)arg8.data_d;
    double *data10 = (double *)arg10.data_d;
    double *data12 = (double *)arg12.data_d;
    double *data14 = (double *)arg14.data_d;
    double *data16 = (double *)arg16.data_d;

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

      #pragma acc parallel loop independent deviceptr(col_reord,map8,data0,data1,data2,data3,data4,data5,data6,data7,data18,data19,data8,data10,data12,data14,data16)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map8idx;
        int map9idx;
        map8idx = map8[n + set_size1 * 0];
        map9idx = map8[n + set_size1 * 1];


        poisson_op2_openacc(
          &data0[2 * n],
          &data1[1 * n],
          &data2[60 * n],
          &data3[60 * n],
          &data4[60 * n],
          &data5[60 * n],
          &data6[60 * n],
          &data7[60 * n],
          &data8[18 * map8idx],
          &data8[18 * map9idx],
          &data10[1 * map8idx],
          &data10[1 * map9idx],
          &data12[18 * map8idx],
          &data12[18 * map9idx],
          &data14[10 * map8idx],
          &data14[10 * map9idx],
          &data16[100 * map8idx],
          &data16[100 * map9idx],
          &data18[100 * n],
          &data19[100 * n]);
      }

    }
    OP_kernels[23].transfer  += Plan->transfer;
    OP_kernels[23].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[23].time     += wall_t2 - wall_t1;
}
