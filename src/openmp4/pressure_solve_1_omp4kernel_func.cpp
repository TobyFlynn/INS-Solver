//
// auto-generated by op2.py
//

void pressure_solve_1_omp4_kernel(
  int *data0,
  int dat0size,
  bool *data1,
  int dat1size,
  int *map2,
  int map2size,
  double *data2,
  int dat2size,
  double *data4,
  int dat4size,
  double *data6,
  int dat6size,
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
  double *data18,
  int dat18size,
  double *data20,
  int dat20size,
  double *data22,
  int dat22size,
  double *data24,
  int dat24size,
  double *data26,
  int dat26size,
  double *data28,
  int dat28size,
  double *data30,
  int dat30size,
  double *data32,
  int dat32size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size]) \
    map(to: gaussW_g_ompkernel[:7], gFInterp0_g_ompkernel[:105], gFInterp1_g_ompkernel[:105], gFInterp2_g_ompkernel[:105])\
    map(to:col_reord[0:set_size1],map2[0:map2size],data2[0:dat2size],data4[0:dat4size],data6[0:dat6size],data8[0:dat8size],data10[0:dat10size],data12[0:dat12size],data14[0:dat14size],data16[0:dat16size],data18[0:dat18size],data20[0:dat20size],data22[0:dat22size],data24[0:dat24size],data26[0:dat26size],data28[0:dat28size],data30[0:dat30size],data32[0:dat32size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map2idx;
    int map3idx;
    map2idx = map2[n_op + set_size1 * 0];
    map3idx = map2[n_op + set_size1 * 1];

    const double* arg2_vec[] = {
       &data2[105 * map2idx],
       &data2[105 * map3idx]};
    const double* arg4_vec[] = {
       &data4[105 * map2idx],
       &data4[105 * map3idx]};
    const double* arg6_vec[] = {
       &data6[105 * map2idx],
       &data6[105 * map3idx]};
    const double* arg8_vec[] = {
       &data8[105 * map2idx],
       &data8[105 * map3idx]};
    const double* arg10_vec[] = {
       &data10[105 * map2idx],
       &data10[105 * map3idx]};
    const double* arg12_vec[] = {
       &data12[105 * map2idx],
       &data12[105 * map3idx]};
    const double* arg14_vec[] = {
       &data14[105 * map2idx],
       &data14[105 * map3idx]};
    const double* arg16_vec[] = {
       &data16[105 * map2idx],
       &data16[105 * map3idx]};
    const double* arg18_vec[] = {
       &data18[105 * map2idx],
       &data18[105 * map3idx]};
    const double* arg20_vec[] = {
       &data20[21 * map2idx],
       &data20[21 * map3idx]};
    const double* arg22_vec[] = {
       &data22[1 * map2idx],
       &data22[1 * map3idx]};
    const double* arg24_vec[] = {
       &data24[3 * map2idx],
       &data24[3 * map3idx]};
    const double* arg26_vec[] = {
       &data26[21 * map2idx],
       &data26[21 * map3idx]};
    const double* arg28_vec[] = {
       &data28[15 * map2idx],
       &data28[15 * map3idx]};
    const double* arg30_vec[] = {
       &data30[15 * map2idx],
       &data30[15 * map3idx]};
    //variable mapping
    const int *edgeNum = &data0[2*n_op];
    const bool *rev = &data1[1*n_op];
    const double **mD0 = arg2_vec;
    const double **mD1 = arg4_vec;
    const double **mD2 = arg6_vec;
    const double **pD0 = arg8_vec;
    const double **pD1 = arg10_vec;
    const double **pD2 = arg12_vec;
    const double **gVP0 = arg14_vec;
    const double **gVP1 = arg16_vec;
    const double **gVP2 = arg18_vec;
    const double **sJ = arg20_vec;
    const double **h = arg22_vec;
    const double **tau = arg24_vec;
    const double **gRho = arg26_vec;
    const double **rho = arg28_vec;
    const double **u = arg30_vec;
    double *rhsL = &data32[15 * map2idx];
    double *rhsR = &data32[15 * map3idx];

    //inline function
    

    int edgeL = edgeNum[0];
    int edgeR = edgeNum[1];
    bool reverse = *rev;


    const double *mDL, *mDR, *pDL, *pDR, *gVML, *gVMR, *gVPL, *gVPR;
    if(edgeL == 0) {
      mDL  = mD0[0];
      pDL  = pD0[0];
      gVML = gFInterp0_g_ompkernel;
      gVPL = gVP0[0];
    } else if(edgeL == 1) {
      mDL  = mD1[0];
      pDL  = pD1[0];
      gVML = gFInterp1_g_ompkernel;
      gVPL = gVP1[0];
    } else {
      mDL  = mD2[0];
      pDL  = pD2[0];
      gVML = gFInterp2_g_ompkernel;
      gVPL = gVP2[0];
    }

    if(edgeR == 0) {
      mDR  = mD0[1];
      pDR  = pD0[1];
      gVMR = gFInterp0_g_ompkernel;
      gVPR = gVP0[1];
    } else if(edgeR == 1) {
      mDR  = mD1[1];
      pDR  = pD1[1];
      gVMR = gFInterp1_g_ompkernel;
      gVPR = gVP1[1];
    } else {
      mDR  = mD2[1];
      pDR  = pD2[1];
      gVMR = gFInterp2_g_ompkernel;
      gVPR = gVP2[1];
    }

    double op1L[15 * 15];
    double op1R[15 * 15];
    double op2L[15 * 15];
    double op2R[15 * 15];



    for(int i = 0; i < 15; i++) {
      for(int j = 0; j < 15; j++) {
        int c_ind = i * 15 + j;
        op1L[c_ind] = 0.0;
        op1R[c_ind] = 0.0;
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

          op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_ompkernel[k] * sJ[0][factors_indL]
                         * (1.0 / gRho[0][factors_indL]) * mDL[b_ind];
          op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_ompkernel[k] * sJ[1][factors_indR]
                         * (1.0 / gRho[1][factors_indR]) * mDR[b_ind];

          op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_ompkernel[k] * sJ[0][factors_indL]
                         * (1.0 / gRho[1][factors_indRR]) * pDL[b_ind];
          op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_ompkernel[k] * sJ[1][factors_indR]
                         * (1.0 / gRho[0][factors_indLR]) * pDR[b_ind];









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










          op1L[c_ind] += -(1.0 / gRho[0][factors_indL]) * mDL[a_ind] * gaussW_g_ompkernel[k]
                         * sJ[0][factors_indL] * gVML[b_ind];
          op1R[c_ind] += -(1.0 / gRho[1][factors_indR]) * mDR[a_ind] * gaussW_g_ompkernel[k]
                         * sJ[1][factors_indR] * gVMR[b_ind];

          op2L[c_ind] += (1.0 / gRho[0][factors_indL]) * mDL[a_ind] * gaussW_g_ompkernel[k]
                         * sJ[0][factors_indL] * gVPL[b_ind];
          op2R[c_ind] += (1.0 / gRho[1][factors_indR]) * mDR[a_ind] * gaussW_g_ompkernel[k]
                         * sJ[1][factors_indR] * gVPR[b_ind];









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
      tauL[i] = 10 * 0.5 * 5 * 6 * fmax(*(h[0]) / gRho[0][indL], *(h[1]) / gRho[1][indR]);
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
      tauR[i] = 10 * 0.5 * 5 * 6 * fmax(*(h[0]) / gRho[0][indL], *(h[1]) / gRho[1][indR]);
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

          op1L[c_ind] += gVML[a_ind] * gaussW_g_ompkernel[k] * sJ[0][factors_indL]
                         * tauL[k] * gVML[b_ind];
          op1R[c_ind] += gVMR[a_ind] * gaussW_g_ompkernel[k] * sJ[1][factors_indR]
                         * tauR[k] * gVMR[b_ind];

          op2L[c_ind] += -gVML[a_ind] * gaussW_g_ompkernel[k] * sJ[0][factors_indL]
                         * tauL[k] * gVPL[b_ind];
          op2R[c_ind] += -gVMR[a_ind] * gaussW_g_ompkernel[k] * sJ[1][factors_indR]
                         * tauR[k] * gVPR[b_ind];









        }
      }
    }


    for(int i = 0; i < 15; i++) {
      for(int j = 0; j < 15; j++) {
        int op_ind = i * 15 + j;
        rhsL[i] += op1L[op_ind] * u[0][j] + op2L[op_ind] * u[1][j];
        rhsR[i] += op1R[op_ind] * u[1][j] + op2R[op_ind] * u[0][j];
      }
    }
    //end inline func
  }

}
