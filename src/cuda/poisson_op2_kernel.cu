//
// auto-generated by op2.py
//

//user function
__device__ void poisson_op2_gpu( const int *edgeNum, const bool *rev,
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
    gVML = gFInterp0_g_cuda;
    gVPL = gVP0L;
  } else if(edgeL == 1) {
    mDL  = mD1L;
    pDL  = pD1L;
    gVML = gFInterp1_g_cuda;
    gVPL = gVP1L;
  } else {
    mDL  = mD2L;
    pDL  = pD2L;
    gVML = gFInterp2_g_cuda;
    gVPL = gVP2L;
  }

  if(edgeR == 0) {
    mDR  = mD0R;
    pDR  = pD0R;
    gVMR = gFInterp0_g_cuda;
    gVPR = gVP0R;
  } else if(edgeR == 1) {
    mDR  = mD1R;
    pDR  = pD1R;
    gVMR = gFInterp1_g_cuda;
    gVPR = gVP1R;
  } else {
    mDR  = mD2R;
    pDR  = pD2R;
    gVMR = gFInterp2_g_cuda;
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

        op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * gFactorL[factors_indL] * mDL[b_ind];
        op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
                       * gFactorR[factors_indR] * mDR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * gFactorR[factors_indRR] * pDL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
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

        op1L[c_ind] += -factorL[i] * mDL[a_ind] * gaussW_g_cuda[k]
                       * sJL[factors_indL] * gVML[b_ind];
        op1R[c_ind] += -factorR[i] * mDR[a_ind] * gaussW_g_cuda[k]
                       * sJR[factors_indR] * gVMR[b_ind];

        op2L[c_ind] += factorL[i] * mDL[a_ind] * gaussW_g_cuda[k]
                       * sJL[factors_indL] * gVPL[b_ind];
        op2R[c_ind] += factorR[i] * mDR[a_ind] * gaussW_g_cuda[k]
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

        op1L[c_ind] += gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * tauL[k] * gVML[b_ind];
        op1R[c_ind] += gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
                       * tauR[k] * gVMR[b_ind];

        op2L[c_ind] += -gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * tauL[k] * gVPL[b_ind];
        op2R[c_ind] += -gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
                       * tauR[k] * gVPR[b_ind];
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_op2(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  const double *__restrict ind_arg6,
  const double *__restrict ind_arg7,
  const double *__restrict ind_arg8,
  const double *__restrict ind_arg9,
  const double *__restrict ind_arg10,
  const double *__restrict ind_arg11,
  const double *__restrict ind_arg12,
  double *__restrict ind_arg13,
  const int *__restrict opDat2Map,
  const int *__restrict arg0,
  const bool *__restrict arg1,
  double *arg30,
  double *arg31,
  int start,
  int end,
  int   set_size) {
  double arg28_l[225];
  double arg29_l[225];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg28_l[225];
    for ( int d=0; d<225; d++ ){
      arg28_l[d] = ZERO_double;
    }
    double arg29_l[225];
    for ( int d=0; d<225; d++ ){
      arg29_l[d] = ZERO_double;
    }
    int map2idx;
    int map3idx;
    map2idx = opDat2Map[n + set_size * 0];
    map3idx = opDat2Map[n + set_size * 1];

    //user-supplied kernel call
    poisson_op2_gpu(arg0+n*2,
                arg1+n*1,
                ind_arg0+map2idx*105,
                ind_arg0+map3idx*105,
                ind_arg1+map2idx*105,
                ind_arg1+map3idx*105,
                ind_arg2+map2idx*105,
                ind_arg2+map3idx*105,
                ind_arg3+map2idx*105,
                ind_arg3+map3idx*105,
                ind_arg4+map2idx*105,
                ind_arg4+map3idx*105,
                ind_arg5+map2idx*105,
                ind_arg5+map3idx*105,
                ind_arg6+map2idx*105,
                ind_arg6+map3idx*105,
                ind_arg7+map2idx*105,
                ind_arg7+map3idx*105,
                ind_arg8+map2idx*105,
                ind_arg8+map3idx*105,
                ind_arg9+map2idx*21,
                ind_arg9+map3idx*21,
                ind_arg10+map2idx*1,
                ind_arg10+map3idx*1,
                ind_arg11+map2idx*21,
                ind_arg11+map3idx*21,
                ind_arg12+map2idx*15,
                ind_arg12+map3idx*15,
                arg28_l,
                arg29_l,
                arg30+n*225,
                arg31+n*225);
    atomicAdd(&ind_arg13[0+map2idx*225],arg28_l[0]);
    atomicAdd(&ind_arg13[1+map2idx*225],arg28_l[1]);
    atomicAdd(&ind_arg13[2+map2idx*225],arg28_l[2]);
    atomicAdd(&ind_arg13[3+map2idx*225],arg28_l[3]);
    atomicAdd(&ind_arg13[4+map2idx*225],arg28_l[4]);
    atomicAdd(&ind_arg13[5+map2idx*225],arg28_l[5]);
    atomicAdd(&ind_arg13[6+map2idx*225],arg28_l[6]);
    atomicAdd(&ind_arg13[7+map2idx*225],arg28_l[7]);
    atomicAdd(&ind_arg13[8+map2idx*225],arg28_l[8]);
    atomicAdd(&ind_arg13[9+map2idx*225],arg28_l[9]);
    atomicAdd(&ind_arg13[10+map2idx*225],arg28_l[10]);
    atomicAdd(&ind_arg13[11+map2idx*225],arg28_l[11]);
    atomicAdd(&ind_arg13[12+map2idx*225],arg28_l[12]);
    atomicAdd(&ind_arg13[13+map2idx*225],arg28_l[13]);
    atomicAdd(&ind_arg13[14+map2idx*225],arg28_l[14]);
    atomicAdd(&ind_arg13[15+map2idx*225],arg28_l[15]);
    atomicAdd(&ind_arg13[16+map2idx*225],arg28_l[16]);
    atomicAdd(&ind_arg13[17+map2idx*225],arg28_l[17]);
    atomicAdd(&ind_arg13[18+map2idx*225],arg28_l[18]);
    atomicAdd(&ind_arg13[19+map2idx*225],arg28_l[19]);
    atomicAdd(&ind_arg13[20+map2idx*225],arg28_l[20]);
    atomicAdd(&ind_arg13[21+map2idx*225],arg28_l[21]);
    atomicAdd(&ind_arg13[22+map2idx*225],arg28_l[22]);
    atomicAdd(&ind_arg13[23+map2idx*225],arg28_l[23]);
    atomicAdd(&ind_arg13[24+map2idx*225],arg28_l[24]);
    atomicAdd(&ind_arg13[25+map2idx*225],arg28_l[25]);
    atomicAdd(&ind_arg13[26+map2idx*225],arg28_l[26]);
    atomicAdd(&ind_arg13[27+map2idx*225],arg28_l[27]);
    atomicAdd(&ind_arg13[28+map2idx*225],arg28_l[28]);
    atomicAdd(&ind_arg13[29+map2idx*225],arg28_l[29]);
    atomicAdd(&ind_arg13[30+map2idx*225],arg28_l[30]);
    atomicAdd(&ind_arg13[31+map2idx*225],arg28_l[31]);
    atomicAdd(&ind_arg13[32+map2idx*225],arg28_l[32]);
    atomicAdd(&ind_arg13[33+map2idx*225],arg28_l[33]);
    atomicAdd(&ind_arg13[34+map2idx*225],arg28_l[34]);
    atomicAdd(&ind_arg13[35+map2idx*225],arg28_l[35]);
    atomicAdd(&ind_arg13[36+map2idx*225],arg28_l[36]);
    atomicAdd(&ind_arg13[37+map2idx*225],arg28_l[37]);
    atomicAdd(&ind_arg13[38+map2idx*225],arg28_l[38]);
    atomicAdd(&ind_arg13[39+map2idx*225],arg28_l[39]);
    atomicAdd(&ind_arg13[40+map2idx*225],arg28_l[40]);
    atomicAdd(&ind_arg13[41+map2idx*225],arg28_l[41]);
    atomicAdd(&ind_arg13[42+map2idx*225],arg28_l[42]);
    atomicAdd(&ind_arg13[43+map2idx*225],arg28_l[43]);
    atomicAdd(&ind_arg13[44+map2idx*225],arg28_l[44]);
    atomicAdd(&ind_arg13[45+map2idx*225],arg28_l[45]);
    atomicAdd(&ind_arg13[46+map2idx*225],arg28_l[46]);
    atomicAdd(&ind_arg13[47+map2idx*225],arg28_l[47]);
    atomicAdd(&ind_arg13[48+map2idx*225],arg28_l[48]);
    atomicAdd(&ind_arg13[49+map2idx*225],arg28_l[49]);
    atomicAdd(&ind_arg13[50+map2idx*225],arg28_l[50]);
    atomicAdd(&ind_arg13[51+map2idx*225],arg28_l[51]);
    atomicAdd(&ind_arg13[52+map2idx*225],arg28_l[52]);
    atomicAdd(&ind_arg13[53+map2idx*225],arg28_l[53]);
    atomicAdd(&ind_arg13[54+map2idx*225],arg28_l[54]);
    atomicAdd(&ind_arg13[55+map2idx*225],arg28_l[55]);
    atomicAdd(&ind_arg13[56+map2idx*225],arg28_l[56]);
    atomicAdd(&ind_arg13[57+map2idx*225],arg28_l[57]);
    atomicAdd(&ind_arg13[58+map2idx*225],arg28_l[58]);
    atomicAdd(&ind_arg13[59+map2idx*225],arg28_l[59]);
    atomicAdd(&ind_arg13[60+map2idx*225],arg28_l[60]);
    atomicAdd(&ind_arg13[61+map2idx*225],arg28_l[61]);
    atomicAdd(&ind_arg13[62+map2idx*225],arg28_l[62]);
    atomicAdd(&ind_arg13[63+map2idx*225],arg28_l[63]);
    atomicAdd(&ind_arg13[64+map2idx*225],arg28_l[64]);
    atomicAdd(&ind_arg13[65+map2idx*225],arg28_l[65]);
    atomicAdd(&ind_arg13[66+map2idx*225],arg28_l[66]);
    atomicAdd(&ind_arg13[67+map2idx*225],arg28_l[67]);
    atomicAdd(&ind_arg13[68+map2idx*225],arg28_l[68]);
    atomicAdd(&ind_arg13[69+map2idx*225],arg28_l[69]);
    atomicAdd(&ind_arg13[70+map2idx*225],arg28_l[70]);
    atomicAdd(&ind_arg13[71+map2idx*225],arg28_l[71]);
    atomicAdd(&ind_arg13[72+map2idx*225],arg28_l[72]);
    atomicAdd(&ind_arg13[73+map2idx*225],arg28_l[73]);
    atomicAdd(&ind_arg13[74+map2idx*225],arg28_l[74]);
    atomicAdd(&ind_arg13[75+map2idx*225],arg28_l[75]);
    atomicAdd(&ind_arg13[76+map2idx*225],arg28_l[76]);
    atomicAdd(&ind_arg13[77+map2idx*225],arg28_l[77]);
    atomicAdd(&ind_arg13[78+map2idx*225],arg28_l[78]);
    atomicAdd(&ind_arg13[79+map2idx*225],arg28_l[79]);
    atomicAdd(&ind_arg13[80+map2idx*225],arg28_l[80]);
    atomicAdd(&ind_arg13[81+map2idx*225],arg28_l[81]);
    atomicAdd(&ind_arg13[82+map2idx*225],arg28_l[82]);
    atomicAdd(&ind_arg13[83+map2idx*225],arg28_l[83]);
    atomicAdd(&ind_arg13[84+map2idx*225],arg28_l[84]);
    atomicAdd(&ind_arg13[85+map2idx*225],arg28_l[85]);
    atomicAdd(&ind_arg13[86+map2idx*225],arg28_l[86]);
    atomicAdd(&ind_arg13[87+map2idx*225],arg28_l[87]);
    atomicAdd(&ind_arg13[88+map2idx*225],arg28_l[88]);
    atomicAdd(&ind_arg13[89+map2idx*225],arg28_l[89]);
    atomicAdd(&ind_arg13[90+map2idx*225],arg28_l[90]);
    atomicAdd(&ind_arg13[91+map2idx*225],arg28_l[91]);
    atomicAdd(&ind_arg13[92+map2idx*225],arg28_l[92]);
    atomicAdd(&ind_arg13[93+map2idx*225],arg28_l[93]);
    atomicAdd(&ind_arg13[94+map2idx*225],arg28_l[94]);
    atomicAdd(&ind_arg13[95+map2idx*225],arg28_l[95]);
    atomicAdd(&ind_arg13[96+map2idx*225],arg28_l[96]);
    atomicAdd(&ind_arg13[97+map2idx*225],arg28_l[97]);
    atomicAdd(&ind_arg13[98+map2idx*225],arg28_l[98]);
    atomicAdd(&ind_arg13[99+map2idx*225],arg28_l[99]);
    atomicAdd(&ind_arg13[100+map2idx*225],arg28_l[100]);
    atomicAdd(&ind_arg13[101+map2idx*225],arg28_l[101]);
    atomicAdd(&ind_arg13[102+map2idx*225],arg28_l[102]);
    atomicAdd(&ind_arg13[103+map2idx*225],arg28_l[103]);
    atomicAdd(&ind_arg13[104+map2idx*225],arg28_l[104]);
    atomicAdd(&ind_arg13[105+map2idx*225],arg28_l[105]);
    atomicAdd(&ind_arg13[106+map2idx*225],arg28_l[106]);
    atomicAdd(&ind_arg13[107+map2idx*225],arg28_l[107]);
    atomicAdd(&ind_arg13[108+map2idx*225],arg28_l[108]);
    atomicAdd(&ind_arg13[109+map2idx*225],arg28_l[109]);
    atomicAdd(&ind_arg13[110+map2idx*225],arg28_l[110]);
    atomicAdd(&ind_arg13[111+map2idx*225],arg28_l[111]);
    atomicAdd(&ind_arg13[112+map2idx*225],arg28_l[112]);
    atomicAdd(&ind_arg13[113+map2idx*225],arg28_l[113]);
    atomicAdd(&ind_arg13[114+map2idx*225],arg28_l[114]);
    atomicAdd(&ind_arg13[115+map2idx*225],arg28_l[115]);
    atomicAdd(&ind_arg13[116+map2idx*225],arg28_l[116]);
    atomicAdd(&ind_arg13[117+map2idx*225],arg28_l[117]);
    atomicAdd(&ind_arg13[118+map2idx*225],arg28_l[118]);
    atomicAdd(&ind_arg13[119+map2idx*225],arg28_l[119]);
    atomicAdd(&ind_arg13[120+map2idx*225],arg28_l[120]);
    atomicAdd(&ind_arg13[121+map2idx*225],arg28_l[121]);
    atomicAdd(&ind_arg13[122+map2idx*225],arg28_l[122]);
    atomicAdd(&ind_arg13[123+map2idx*225],arg28_l[123]);
    atomicAdd(&ind_arg13[124+map2idx*225],arg28_l[124]);
    atomicAdd(&ind_arg13[125+map2idx*225],arg28_l[125]);
    atomicAdd(&ind_arg13[126+map2idx*225],arg28_l[126]);
    atomicAdd(&ind_arg13[127+map2idx*225],arg28_l[127]);
    atomicAdd(&ind_arg13[128+map2idx*225],arg28_l[128]);
    atomicAdd(&ind_arg13[129+map2idx*225],arg28_l[129]);
    atomicAdd(&ind_arg13[130+map2idx*225],arg28_l[130]);
    atomicAdd(&ind_arg13[131+map2idx*225],arg28_l[131]);
    atomicAdd(&ind_arg13[132+map2idx*225],arg28_l[132]);
    atomicAdd(&ind_arg13[133+map2idx*225],arg28_l[133]);
    atomicAdd(&ind_arg13[134+map2idx*225],arg28_l[134]);
    atomicAdd(&ind_arg13[135+map2idx*225],arg28_l[135]);
    atomicAdd(&ind_arg13[136+map2idx*225],arg28_l[136]);
    atomicAdd(&ind_arg13[137+map2idx*225],arg28_l[137]);
    atomicAdd(&ind_arg13[138+map2idx*225],arg28_l[138]);
    atomicAdd(&ind_arg13[139+map2idx*225],arg28_l[139]);
    atomicAdd(&ind_arg13[140+map2idx*225],arg28_l[140]);
    atomicAdd(&ind_arg13[141+map2idx*225],arg28_l[141]);
    atomicAdd(&ind_arg13[142+map2idx*225],arg28_l[142]);
    atomicAdd(&ind_arg13[143+map2idx*225],arg28_l[143]);
    atomicAdd(&ind_arg13[144+map2idx*225],arg28_l[144]);
    atomicAdd(&ind_arg13[145+map2idx*225],arg28_l[145]);
    atomicAdd(&ind_arg13[146+map2idx*225],arg28_l[146]);
    atomicAdd(&ind_arg13[147+map2idx*225],arg28_l[147]);
    atomicAdd(&ind_arg13[148+map2idx*225],arg28_l[148]);
    atomicAdd(&ind_arg13[149+map2idx*225],arg28_l[149]);
    atomicAdd(&ind_arg13[150+map2idx*225],arg28_l[150]);
    atomicAdd(&ind_arg13[151+map2idx*225],arg28_l[151]);
    atomicAdd(&ind_arg13[152+map2idx*225],arg28_l[152]);
    atomicAdd(&ind_arg13[153+map2idx*225],arg28_l[153]);
    atomicAdd(&ind_arg13[154+map2idx*225],arg28_l[154]);
    atomicAdd(&ind_arg13[155+map2idx*225],arg28_l[155]);
    atomicAdd(&ind_arg13[156+map2idx*225],arg28_l[156]);
    atomicAdd(&ind_arg13[157+map2idx*225],arg28_l[157]);
    atomicAdd(&ind_arg13[158+map2idx*225],arg28_l[158]);
    atomicAdd(&ind_arg13[159+map2idx*225],arg28_l[159]);
    atomicAdd(&ind_arg13[160+map2idx*225],arg28_l[160]);
    atomicAdd(&ind_arg13[161+map2idx*225],arg28_l[161]);
    atomicAdd(&ind_arg13[162+map2idx*225],arg28_l[162]);
    atomicAdd(&ind_arg13[163+map2idx*225],arg28_l[163]);
    atomicAdd(&ind_arg13[164+map2idx*225],arg28_l[164]);
    atomicAdd(&ind_arg13[165+map2idx*225],arg28_l[165]);
    atomicAdd(&ind_arg13[166+map2idx*225],arg28_l[166]);
    atomicAdd(&ind_arg13[167+map2idx*225],arg28_l[167]);
    atomicAdd(&ind_arg13[168+map2idx*225],arg28_l[168]);
    atomicAdd(&ind_arg13[169+map2idx*225],arg28_l[169]);
    atomicAdd(&ind_arg13[170+map2idx*225],arg28_l[170]);
    atomicAdd(&ind_arg13[171+map2idx*225],arg28_l[171]);
    atomicAdd(&ind_arg13[172+map2idx*225],arg28_l[172]);
    atomicAdd(&ind_arg13[173+map2idx*225],arg28_l[173]);
    atomicAdd(&ind_arg13[174+map2idx*225],arg28_l[174]);
    atomicAdd(&ind_arg13[175+map2idx*225],arg28_l[175]);
    atomicAdd(&ind_arg13[176+map2idx*225],arg28_l[176]);
    atomicAdd(&ind_arg13[177+map2idx*225],arg28_l[177]);
    atomicAdd(&ind_arg13[178+map2idx*225],arg28_l[178]);
    atomicAdd(&ind_arg13[179+map2idx*225],arg28_l[179]);
    atomicAdd(&ind_arg13[180+map2idx*225],arg28_l[180]);
    atomicAdd(&ind_arg13[181+map2idx*225],arg28_l[181]);
    atomicAdd(&ind_arg13[182+map2idx*225],arg28_l[182]);
    atomicAdd(&ind_arg13[183+map2idx*225],arg28_l[183]);
    atomicAdd(&ind_arg13[184+map2idx*225],arg28_l[184]);
    atomicAdd(&ind_arg13[185+map2idx*225],arg28_l[185]);
    atomicAdd(&ind_arg13[186+map2idx*225],arg28_l[186]);
    atomicAdd(&ind_arg13[187+map2idx*225],arg28_l[187]);
    atomicAdd(&ind_arg13[188+map2idx*225],arg28_l[188]);
    atomicAdd(&ind_arg13[189+map2idx*225],arg28_l[189]);
    atomicAdd(&ind_arg13[190+map2idx*225],arg28_l[190]);
    atomicAdd(&ind_arg13[191+map2idx*225],arg28_l[191]);
    atomicAdd(&ind_arg13[192+map2idx*225],arg28_l[192]);
    atomicAdd(&ind_arg13[193+map2idx*225],arg28_l[193]);
    atomicAdd(&ind_arg13[194+map2idx*225],arg28_l[194]);
    atomicAdd(&ind_arg13[195+map2idx*225],arg28_l[195]);
    atomicAdd(&ind_arg13[196+map2idx*225],arg28_l[196]);
    atomicAdd(&ind_arg13[197+map2idx*225],arg28_l[197]);
    atomicAdd(&ind_arg13[198+map2idx*225],arg28_l[198]);
    atomicAdd(&ind_arg13[199+map2idx*225],arg28_l[199]);
    atomicAdd(&ind_arg13[200+map2idx*225],arg28_l[200]);
    atomicAdd(&ind_arg13[201+map2idx*225],arg28_l[201]);
    atomicAdd(&ind_arg13[202+map2idx*225],arg28_l[202]);
    atomicAdd(&ind_arg13[203+map2idx*225],arg28_l[203]);
    atomicAdd(&ind_arg13[204+map2idx*225],arg28_l[204]);
    atomicAdd(&ind_arg13[205+map2idx*225],arg28_l[205]);
    atomicAdd(&ind_arg13[206+map2idx*225],arg28_l[206]);
    atomicAdd(&ind_arg13[207+map2idx*225],arg28_l[207]);
    atomicAdd(&ind_arg13[208+map2idx*225],arg28_l[208]);
    atomicAdd(&ind_arg13[209+map2idx*225],arg28_l[209]);
    atomicAdd(&ind_arg13[210+map2idx*225],arg28_l[210]);
    atomicAdd(&ind_arg13[211+map2idx*225],arg28_l[211]);
    atomicAdd(&ind_arg13[212+map2idx*225],arg28_l[212]);
    atomicAdd(&ind_arg13[213+map2idx*225],arg28_l[213]);
    atomicAdd(&ind_arg13[214+map2idx*225],arg28_l[214]);
    atomicAdd(&ind_arg13[215+map2idx*225],arg28_l[215]);
    atomicAdd(&ind_arg13[216+map2idx*225],arg28_l[216]);
    atomicAdd(&ind_arg13[217+map2idx*225],arg28_l[217]);
    atomicAdd(&ind_arg13[218+map2idx*225],arg28_l[218]);
    atomicAdd(&ind_arg13[219+map2idx*225],arg28_l[219]);
    atomicAdd(&ind_arg13[220+map2idx*225],arg28_l[220]);
    atomicAdd(&ind_arg13[221+map2idx*225],arg28_l[221]);
    atomicAdd(&ind_arg13[222+map2idx*225],arg28_l[222]);
    atomicAdd(&ind_arg13[223+map2idx*225],arg28_l[223]);
    atomicAdd(&ind_arg13[224+map2idx*225],arg28_l[224]);
    atomicAdd(&ind_arg13[0+map3idx*225],arg29_l[0]);
    atomicAdd(&ind_arg13[1+map3idx*225],arg29_l[1]);
    atomicAdd(&ind_arg13[2+map3idx*225],arg29_l[2]);
    atomicAdd(&ind_arg13[3+map3idx*225],arg29_l[3]);
    atomicAdd(&ind_arg13[4+map3idx*225],arg29_l[4]);
    atomicAdd(&ind_arg13[5+map3idx*225],arg29_l[5]);
    atomicAdd(&ind_arg13[6+map3idx*225],arg29_l[6]);
    atomicAdd(&ind_arg13[7+map3idx*225],arg29_l[7]);
    atomicAdd(&ind_arg13[8+map3idx*225],arg29_l[8]);
    atomicAdd(&ind_arg13[9+map3idx*225],arg29_l[9]);
    atomicAdd(&ind_arg13[10+map3idx*225],arg29_l[10]);
    atomicAdd(&ind_arg13[11+map3idx*225],arg29_l[11]);
    atomicAdd(&ind_arg13[12+map3idx*225],arg29_l[12]);
    atomicAdd(&ind_arg13[13+map3idx*225],arg29_l[13]);
    atomicAdd(&ind_arg13[14+map3idx*225],arg29_l[14]);
    atomicAdd(&ind_arg13[15+map3idx*225],arg29_l[15]);
    atomicAdd(&ind_arg13[16+map3idx*225],arg29_l[16]);
    atomicAdd(&ind_arg13[17+map3idx*225],arg29_l[17]);
    atomicAdd(&ind_arg13[18+map3idx*225],arg29_l[18]);
    atomicAdd(&ind_arg13[19+map3idx*225],arg29_l[19]);
    atomicAdd(&ind_arg13[20+map3idx*225],arg29_l[20]);
    atomicAdd(&ind_arg13[21+map3idx*225],arg29_l[21]);
    atomicAdd(&ind_arg13[22+map3idx*225],arg29_l[22]);
    atomicAdd(&ind_arg13[23+map3idx*225],arg29_l[23]);
    atomicAdd(&ind_arg13[24+map3idx*225],arg29_l[24]);
    atomicAdd(&ind_arg13[25+map3idx*225],arg29_l[25]);
    atomicAdd(&ind_arg13[26+map3idx*225],arg29_l[26]);
    atomicAdd(&ind_arg13[27+map3idx*225],arg29_l[27]);
    atomicAdd(&ind_arg13[28+map3idx*225],arg29_l[28]);
    atomicAdd(&ind_arg13[29+map3idx*225],arg29_l[29]);
    atomicAdd(&ind_arg13[30+map3idx*225],arg29_l[30]);
    atomicAdd(&ind_arg13[31+map3idx*225],arg29_l[31]);
    atomicAdd(&ind_arg13[32+map3idx*225],arg29_l[32]);
    atomicAdd(&ind_arg13[33+map3idx*225],arg29_l[33]);
    atomicAdd(&ind_arg13[34+map3idx*225],arg29_l[34]);
    atomicAdd(&ind_arg13[35+map3idx*225],arg29_l[35]);
    atomicAdd(&ind_arg13[36+map3idx*225],arg29_l[36]);
    atomicAdd(&ind_arg13[37+map3idx*225],arg29_l[37]);
    atomicAdd(&ind_arg13[38+map3idx*225],arg29_l[38]);
    atomicAdd(&ind_arg13[39+map3idx*225],arg29_l[39]);
    atomicAdd(&ind_arg13[40+map3idx*225],arg29_l[40]);
    atomicAdd(&ind_arg13[41+map3idx*225],arg29_l[41]);
    atomicAdd(&ind_arg13[42+map3idx*225],arg29_l[42]);
    atomicAdd(&ind_arg13[43+map3idx*225],arg29_l[43]);
    atomicAdd(&ind_arg13[44+map3idx*225],arg29_l[44]);
    atomicAdd(&ind_arg13[45+map3idx*225],arg29_l[45]);
    atomicAdd(&ind_arg13[46+map3idx*225],arg29_l[46]);
    atomicAdd(&ind_arg13[47+map3idx*225],arg29_l[47]);
    atomicAdd(&ind_arg13[48+map3idx*225],arg29_l[48]);
    atomicAdd(&ind_arg13[49+map3idx*225],arg29_l[49]);
    atomicAdd(&ind_arg13[50+map3idx*225],arg29_l[50]);
    atomicAdd(&ind_arg13[51+map3idx*225],arg29_l[51]);
    atomicAdd(&ind_arg13[52+map3idx*225],arg29_l[52]);
    atomicAdd(&ind_arg13[53+map3idx*225],arg29_l[53]);
    atomicAdd(&ind_arg13[54+map3idx*225],arg29_l[54]);
    atomicAdd(&ind_arg13[55+map3idx*225],arg29_l[55]);
    atomicAdd(&ind_arg13[56+map3idx*225],arg29_l[56]);
    atomicAdd(&ind_arg13[57+map3idx*225],arg29_l[57]);
    atomicAdd(&ind_arg13[58+map3idx*225],arg29_l[58]);
    atomicAdd(&ind_arg13[59+map3idx*225],arg29_l[59]);
    atomicAdd(&ind_arg13[60+map3idx*225],arg29_l[60]);
    atomicAdd(&ind_arg13[61+map3idx*225],arg29_l[61]);
    atomicAdd(&ind_arg13[62+map3idx*225],arg29_l[62]);
    atomicAdd(&ind_arg13[63+map3idx*225],arg29_l[63]);
    atomicAdd(&ind_arg13[64+map3idx*225],arg29_l[64]);
    atomicAdd(&ind_arg13[65+map3idx*225],arg29_l[65]);
    atomicAdd(&ind_arg13[66+map3idx*225],arg29_l[66]);
    atomicAdd(&ind_arg13[67+map3idx*225],arg29_l[67]);
    atomicAdd(&ind_arg13[68+map3idx*225],arg29_l[68]);
    atomicAdd(&ind_arg13[69+map3idx*225],arg29_l[69]);
    atomicAdd(&ind_arg13[70+map3idx*225],arg29_l[70]);
    atomicAdd(&ind_arg13[71+map3idx*225],arg29_l[71]);
    atomicAdd(&ind_arg13[72+map3idx*225],arg29_l[72]);
    atomicAdd(&ind_arg13[73+map3idx*225],arg29_l[73]);
    atomicAdd(&ind_arg13[74+map3idx*225],arg29_l[74]);
    atomicAdd(&ind_arg13[75+map3idx*225],arg29_l[75]);
    atomicAdd(&ind_arg13[76+map3idx*225],arg29_l[76]);
    atomicAdd(&ind_arg13[77+map3idx*225],arg29_l[77]);
    atomicAdd(&ind_arg13[78+map3idx*225],arg29_l[78]);
    atomicAdd(&ind_arg13[79+map3idx*225],arg29_l[79]);
    atomicAdd(&ind_arg13[80+map3idx*225],arg29_l[80]);
    atomicAdd(&ind_arg13[81+map3idx*225],arg29_l[81]);
    atomicAdd(&ind_arg13[82+map3idx*225],arg29_l[82]);
    atomicAdd(&ind_arg13[83+map3idx*225],arg29_l[83]);
    atomicAdd(&ind_arg13[84+map3idx*225],arg29_l[84]);
    atomicAdd(&ind_arg13[85+map3idx*225],arg29_l[85]);
    atomicAdd(&ind_arg13[86+map3idx*225],arg29_l[86]);
    atomicAdd(&ind_arg13[87+map3idx*225],arg29_l[87]);
    atomicAdd(&ind_arg13[88+map3idx*225],arg29_l[88]);
    atomicAdd(&ind_arg13[89+map3idx*225],arg29_l[89]);
    atomicAdd(&ind_arg13[90+map3idx*225],arg29_l[90]);
    atomicAdd(&ind_arg13[91+map3idx*225],arg29_l[91]);
    atomicAdd(&ind_arg13[92+map3idx*225],arg29_l[92]);
    atomicAdd(&ind_arg13[93+map3idx*225],arg29_l[93]);
    atomicAdd(&ind_arg13[94+map3idx*225],arg29_l[94]);
    atomicAdd(&ind_arg13[95+map3idx*225],arg29_l[95]);
    atomicAdd(&ind_arg13[96+map3idx*225],arg29_l[96]);
    atomicAdd(&ind_arg13[97+map3idx*225],arg29_l[97]);
    atomicAdd(&ind_arg13[98+map3idx*225],arg29_l[98]);
    atomicAdd(&ind_arg13[99+map3idx*225],arg29_l[99]);
    atomicAdd(&ind_arg13[100+map3idx*225],arg29_l[100]);
    atomicAdd(&ind_arg13[101+map3idx*225],arg29_l[101]);
    atomicAdd(&ind_arg13[102+map3idx*225],arg29_l[102]);
    atomicAdd(&ind_arg13[103+map3idx*225],arg29_l[103]);
    atomicAdd(&ind_arg13[104+map3idx*225],arg29_l[104]);
    atomicAdd(&ind_arg13[105+map3idx*225],arg29_l[105]);
    atomicAdd(&ind_arg13[106+map3idx*225],arg29_l[106]);
    atomicAdd(&ind_arg13[107+map3idx*225],arg29_l[107]);
    atomicAdd(&ind_arg13[108+map3idx*225],arg29_l[108]);
    atomicAdd(&ind_arg13[109+map3idx*225],arg29_l[109]);
    atomicAdd(&ind_arg13[110+map3idx*225],arg29_l[110]);
    atomicAdd(&ind_arg13[111+map3idx*225],arg29_l[111]);
    atomicAdd(&ind_arg13[112+map3idx*225],arg29_l[112]);
    atomicAdd(&ind_arg13[113+map3idx*225],arg29_l[113]);
    atomicAdd(&ind_arg13[114+map3idx*225],arg29_l[114]);
    atomicAdd(&ind_arg13[115+map3idx*225],arg29_l[115]);
    atomicAdd(&ind_arg13[116+map3idx*225],arg29_l[116]);
    atomicAdd(&ind_arg13[117+map3idx*225],arg29_l[117]);
    atomicAdd(&ind_arg13[118+map3idx*225],arg29_l[118]);
    atomicAdd(&ind_arg13[119+map3idx*225],arg29_l[119]);
    atomicAdd(&ind_arg13[120+map3idx*225],arg29_l[120]);
    atomicAdd(&ind_arg13[121+map3idx*225],arg29_l[121]);
    atomicAdd(&ind_arg13[122+map3idx*225],arg29_l[122]);
    atomicAdd(&ind_arg13[123+map3idx*225],arg29_l[123]);
    atomicAdd(&ind_arg13[124+map3idx*225],arg29_l[124]);
    atomicAdd(&ind_arg13[125+map3idx*225],arg29_l[125]);
    atomicAdd(&ind_arg13[126+map3idx*225],arg29_l[126]);
    atomicAdd(&ind_arg13[127+map3idx*225],arg29_l[127]);
    atomicAdd(&ind_arg13[128+map3idx*225],arg29_l[128]);
    atomicAdd(&ind_arg13[129+map3idx*225],arg29_l[129]);
    atomicAdd(&ind_arg13[130+map3idx*225],arg29_l[130]);
    atomicAdd(&ind_arg13[131+map3idx*225],arg29_l[131]);
    atomicAdd(&ind_arg13[132+map3idx*225],arg29_l[132]);
    atomicAdd(&ind_arg13[133+map3idx*225],arg29_l[133]);
    atomicAdd(&ind_arg13[134+map3idx*225],arg29_l[134]);
    atomicAdd(&ind_arg13[135+map3idx*225],arg29_l[135]);
    atomicAdd(&ind_arg13[136+map3idx*225],arg29_l[136]);
    atomicAdd(&ind_arg13[137+map3idx*225],arg29_l[137]);
    atomicAdd(&ind_arg13[138+map3idx*225],arg29_l[138]);
    atomicAdd(&ind_arg13[139+map3idx*225],arg29_l[139]);
    atomicAdd(&ind_arg13[140+map3idx*225],arg29_l[140]);
    atomicAdd(&ind_arg13[141+map3idx*225],arg29_l[141]);
    atomicAdd(&ind_arg13[142+map3idx*225],arg29_l[142]);
    atomicAdd(&ind_arg13[143+map3idx*225],arg29_l[143]);
    atomicAdd(&ind_arg13[144+map3idx*225],arg29_l[144]);
    atomicAdd(&ind_arg13[145+map3idx*225],arg29_l[145]);
    atomicAdd(&ind_arg13[146+map3idx*225],arg29_l[146]);
    atomicAdd(&ind_arg13[147+map3idx*225],arg29_l[147]);
    atomicAdd(&ind_arg13[148+map3idx*225],arg29_l[148]);
    atomicAdd(&ind_arg13[149+map3idx*225],arg29_l[149]);
    atomicAdd(&ind_arg13[150+map3idx*225],arg29_l[150]);
    atomicAdd(&ind_arg13[151+map3idx*225],arg29_l[151]);
    atomicAdd(&ind_arg13[152+map3idx*225],arg29_l[152]);
    atomicAdd(&ind_arg13[153+map3idx*225],arg29_l[153]);
    atomicAdd(&ind_arg13[154+map3idx*225],arg29_l[154]);
    atomicAdd(&ind_arg13[155+map3idx*225],arg29_l[155]);
    atomicAdd(&ind_arg13[156+map3idx*225],arg29_l[156]);
    atomicAdd(&ind_arg13[157+map3idx*225],arg29_l[157]);
    atomicAdd(&ind_arg13[158+map3idx*225],arg29_l[158]);
    atomicAdd(&ind_arg13[159+map3idx*225],arg29_l[159]);
    atomicAdd(&ind_arg13[160+map3idx*225],arg29_l[160]);
    atomicAdd(&ind_arg13[161+map3idx*225],arg29_l[161]);
    atomicAdd(&ind_arg13[162+map3idx*225],arg29_l[162]);
    atomicAdd(&ind_arg13[163+map3idx*225],arg29_l[163]);
    atomicAdd(&ind_arg13[164+map3idx*225],arg29_l[164]);
    atomicAdd(&ind_arg13[165+map3idx*225],arg29_l[165]);
    atomicAdd(&ind_arg13[166+map3idx*225],arg29_l[166]);
    atomicAdd(&ind_arg13[167+map3idx*225],arg29_l[167]);
    atomicAdd(&ind_arg13[168+map3idx*225],arg29_l[168]);
    atomicAdd(&ind_arg13[169+map3idx*225],arg29_l[169]);
    atomicAdd(&ind_arg13[170+map3idx*225],arg29_l[170]);
    atomicAdd(&ind_arg13[171+map3idx*225],arg29_l[171]);
    atomicAdd(&ind_arg13[172+map3idx*225],arg29_l[172]);
    atomicAdd(&ind_arg13[173+map3idx*225],arg29_l[173]);
    atomicAdd(&ind_arg13[174+map3idx*225],arg29_l[174]);
    atomicAdd(&ind_arg13[175+map3idx*225],arg29_l[175]);
    atomicAdd(&ind_arg13[176+map3idx*225],arg29_l[176]);
    atomicAdd(&ind_arg13[177+map3idx*225],arg29_l[177]);
    atomicAdd(&ind_arg13[178+map3idx*225],arg29_l[178]);
    atomicAdd(&ind_arg13[179+map3idx*225],arg29_l[179]);
    atomicAdd(&ind_arg13[180+map3idx*225],arg29_l[180]);
    atomicAdd(&ind_arg13[181+map3idx*225],arg29_l[181]);
    atomicAdd(&ind_arg13[182+map3idx*225],arg29_l[182]);
    atomicAdd(&ind_arg13[183+map3idx*225],arg29_l[183]);
    atomicAdd(&ind_arg13[184+map3idx*225],arg29_l[184]);
    atomicAdd(&ind_arg13[185+map3idx*225],arg29_l[185]);
    atomicAdd(&ind_arg13[186+map3idx*225],arg29_l[186]);
    atomicAdd(&ind_arg13[187+map3idx*225],arg29_l[187]);
    atomicAdd(&ind_arg13[188+map3idx*225],arg29_l[188]);
    atomicAdd(&ind_arg13[189+map3idx*225],arg29_l[189]);
    atomicAdd(&ind_arg13[190+map3idx*225],arg29_l[190]);
    atomicAdd(&ind_arg13[191+map3idx*225],arg29_l[191]);
    atomicAdd(&ind_arg13[192+map3idx*225],arg29_l[192]);
    atomicAdd(&ind_arg13[193+map3idx*225],arg29_l[193]);
    atomicAdd(&ind_arg13[194+map3idx*225],arg29_l[194]);
    atomicAdd(&ind_arg13[195+map3idx*225],arg29_l[195]);
    atomicAdd(&ind_arg13[196+map3idx*225],arg29_l[196]);
    atomicAdd(&ind_arg13[197+map3idx*225],arg29_l[197]);
    atomicAdd(&ind_arg13[198+map3idx*225],arg29_l[198]);
    atomicAdd(&ind_arg13[199+map3idx*225],arg29_l[199]);
    atomicAdd(&ind_arg13[200+map3idx*225],arg29_l[200]);
    atomicAdd(&ind_arg13[201+map3idx*225],arg29_l[201]);
    atomicAdd(&ind_arg13[202+map3idx*225],arg29_l[202]);
    atomicAdd(&ind_arg13[203+map3idx*225],arg29_l[203]);
    atomicAdd(&ind_arg13[204+map3idx*225],arg29_l[204]);
    atomicAdd(&ind_arg13[205+map3idx*225],arg29_l[205]);
    atomicAdd(&ind_arg13[206+map3idx*225],arg29_l[206]);
    atomicAdd(&ind_arg13[207+map3idx*225],arg29_l[207]);
    atomicAdd(&ind_arg13[208+map3idx*225],arg29_l[208]);
    atomicAdd(&ind_arg13[209+map3idx*225],arg29_l[209]);
    atomicAdd(&ind_arg13[210+map3idx*225],arg29_l[210]);
    atomicAdd(&ind_arg13[211+map3idx*225],arg29_l[211]);
    atomicAdd(&ind_arg13[212+map3idx*225],arg29_l[212]);
    atomicAdd(&ind_arg13[213+map3idx*225],arg29_l[213]);
    atomicAdd(&ind_arg13[214+map3idx*225],arg29_l[214]);
    atomicAdd(&ind_arg13[215+map3idx*225],arg29_l[215]);
    atomicAdd(&ind_arg13[216+map3idx*225],arg29_l[216]);
    atomicAdd(&ind_arg13[217+map3idx*225],arg29_l[217]);
    atomicAdd(&ind_arg13[218+map3idx*225],arg29_l[218]);
    atomicAdd(&ind_arg13[219+map3idx*225],arg29_l[219]);
    atomicAdd(&ind_arg13[220+map3idx*225],arg29_l[220]);
    atomicAdd(&ind_arg13[221+map3idx*225],arg29_l[221]);
    atomicAdd(&ind_arg13[222+map3idx*225],arg29_l[222]);
    atomicAdd(&ind_arg13[223+map3idx*225],arg29_l[223]);
    atomicAdd(&ind_arg13[224+map3idx*225],arg29_l[224]);
  }
}


//host stub function
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
  op_timing_realloc(19);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[19].name      = name;
  OP_kernels[19].count    += 1;


  int    ninds   = 14;
  int    inds[32] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_op2\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_19
      int nthread = OP_BLOCK_SIZE_19;
    #else
      int nthread = OP_block_size;
    #endif

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        int nblocks = (end-start-1)/nthread+1;
        op_cuda_poisson_op2<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg4.data_d,
        (double *)arg6.data_d,
        (double *)arg8.data_d,
        (double *)arg10.data_d,
        (double *)arg12.data_d,
        (double *)arg14.data_d,
        (double *)arg16.data_d,
        (double *)arg18.data_d,
        (double *)arg20.data_d,
        (double *)arg22.data_d,
        (double *)arg24.data_d,
        (double *)arg26.data_d,
        (double *)arg28.data_d,
        arg2.map_data_d,
        (int*)arg0.data_d,
        (bool*)arg1.data_d,
        (double*)arg30.data_d,
        (double*)arg31.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[19].time     += wall_t2 - wall_t1;
}
