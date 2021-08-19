//
// auto-generated by op2.py
//

void poisson_op2_omp4_kernel(
  int *data0,
  int dat0size,
  bool *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  double *data5,
  int dat5size,
  double *data6,
  int dat6size,
  double *data7,
  int dat7size,
  int *map8,
  int map8size,
  double *data18,
  int dat18size,
  double *data19,
  int dat19size,
  double *data8,
  int dat8size,
  double *data10,
  int dat10size,
  double *data12,
  int dat12size,
  double *data14,
  int dat14size,
  double *data16,
  int dat16size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data18[0:dat18size],data19[0:dat19size]) \
    map(to: gaussW_g_ompkernel[:6], gFInterp0_g_ompkernel[:60], gFInterp1_g_ompkernel[:60], gFInterp2_g_ompkernel[:60])\
    map(to:col_reord[0:set_size1],map8[0:map8size],data8[0:dat8size],data10[0:dat10size],data12[0:dat12size],data14[0:dat14size],data16[0:dat16size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map8idx;
    int map9idx;
    map8idx = map8[n_op + set_size1 * 0];
    map9idx = map8[n_op + set_size1 * 1];

    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const bool *rev = &data1[1*n_op];
    const double *mDL = &data2[60*n_op];
    const double *mDR = &data3[60*n_op];
    const double *pDL = &data4[60*n_op];
    const double *pDR = &data5[60*n_op];
    const double *gVPL = &data6[60*n_op];
    const double *gVPR = &data7[60*n_op];
    const double *sJL = &data8[18 * map8idx];
    const double *sJR = &data8[18 * map9idx];
    const double *hL = &data10[1 * map8idx];
    const double *hR = &data10[1 * map9idx];
    const double *gFactorL = &data12[18 * map8idx];
    const double *gFactorR = &data12[18 * map9idx];
    const double *factorL = &data14[10 * map8idx];
    const double *factorR = &data14[10 * map9idx];
    double *op1L = &data16[100 * map8idx];
    double *op1R = &data16[100 * map9idx];
    double *op2L = &data18[100*n_op];
    double *op2R = &data19[100*n_op];

    //inline function
    

    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];
    bool reverse = *rev;

    const double *gVML, *gVMR;
    if(edgeL == 0) {
      gVML = gFInterp0_g_ompkernel;
    } else if(edgeL == 1) {
      gVML = gFInterp1_g_ompkernel;
    } else {
      gVML = gFInterp2_g_ompkernel;
    }

    if(edgeR == 0) {
      gVMR = gFInterp0_g_ompkernel;
    } else if(edgeR == 1) {
      gVMR = gFInterp1_g_ompkernel;
    } else {
      gVMR = gFInterp2_g_ompkernel;
    }



    for(int i = 0; i < 10; i++) {
      for(int j = 0; j < 10; j++) {
        int c_ind = i * 10 + j;
        op2L[c_ind] = 0.0;
        op2R[c_ind] = 0.0;
        for(int k = 0; k < 6; k++) {

          int b_ind = k * 10 + j;

          int a_ind = k * 10 + i;

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

          op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_ompkernel[k] * sJL[factors_indL]
                         * gFactorL[factors_indL] * mDL[b_ind];
          op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_ompkernel[k] * sJR[factors_indR]
                         * gFactorR[factors_indR] * mDR[b_ind];

          op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_ompkernel[k] * sJL[factors_indL]
                         * gFactorR[factors_indRR] * pDL[b_ind];
          op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_ompkernel[k] * sJR[factors_indR]
                         * gFactorL[factors_indLR] * pDR[b_ind];
        }
      }
    }



    for(int i = 0; i < 10; i++) {
      for(int j = 0; j < 10; j++) {
        int c_ind = i * 10 + j;
        for(int k = 0; k < 6; k++) {

          int b_ind = k * 10 + j;

          int a_ind = k * 10 + i;

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



















          op1L[c_ind] += -0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g_ompkernel[k]
                         * sJL[factors_indL] * gVML[b_ind];
          op1R[c_ind] += -0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g_ompkernel[k]
                         * sJR[factors_indR] * gVMR[b_ind];

          op2L[c_ind] += 0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g_ompkernel[k]
                         * sJL[factors_indL] * gVPL[b_ind];
          op2R[c_ind] += 0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g_ompkernel[k]
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

          int a_ind = k * 10 + i;

          int factors_indL = edgeL * 6 + k;
          int factors_indR = edgeR * 6 + k;










          op1L[c_ind] += 0.5 * gVML[a_ind] * gaussW_g_ompkernel[k] * sJL[factors_indL]
                         * tauL[k] * gVML[b_ind];
          op1R[c_ind] += 0.5 * gVMR[a_ind] * gaussW_g_ompkernel[k] * sJR[factors_indR]
                         * tauR[k] * gVMR[b_ind];

          op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_ompkernel[k] * sJL[factors_indL]
                         * tauL[k] * gVPL[b_ind];
          op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_ompkernel[k] * sJR[factors_indR]
                         * tauR[k] * gVPR[b_ind];
        }
      }
    }
    //end inline func
  }

}
