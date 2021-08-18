//
// auto-generated by op2.py
//

//user function
__device__ void poisson_op2_gpu( const int *edgeNum, const bool *rev,
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
    gVML = gFInterp0_g_cuda;
  } else if(edgeL == 1) {
    gVML = gFInterp1_g_cuda;
  } else {
    gVML = gFInterp2_g_cuda;
  }

  if(edgeR == 0) {
    gVMR = gFInterp0_g_cuda;
  } else if(edgeR == 1) {
    gVMR = gFInterp1_g_cuda;
  } else {
    gVMR = gFInterp2_g_cuda;
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



















        op1L[c_ind] += -0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g_cuda[k]
                       * sJL[factors_indL] * gVML[b_ind];
        op1R[c_ind] += -0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g_cuda[k]
                       * sJR[factors_indR] * gVMR[b_ind];

        op2L[c_ind] += 0.5 * gFactorL[factors_indL] * mDL[a_ind] * gaussW_g_cuda[k]
                       * sJL[factors_indL] * gVPL[b_ind];
        op2R[c_ind] += 0.5 * gFactorR[factors_indR] * mDR[a_ind] * gaussW_g_cuda[k]
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










        op1L[c_ind] += 0.5 * gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * tauL[k] * gVML[b_ind];
        op1R[c_ind] += 0.5 * gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
                       * tauR[k] * gVMR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * tauL[k] * gVPL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
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
  double *__restrict ind_arg4,
  const int *__restrict opDat8Map,
  const int *__restrict arg0,
  const bool *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  const double *__restrict arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  double *arg18,
  double *arg19,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg16_l[100];
    for ( int d=0; d<100; d++ ){
      arg16_l[d] = ZERO_double;
    }
    double arg17_l[100];
    for ( int d=0; d<100; d++ ){
      arg17_l[d] = ZERO_double;
    }
    int map8idx;
    int map9idx;
    map8idx = opDat8Map[n + set_size * 0];
    map9idx = opDat8Map[n + set_size * 1];

    //user-supplied kernel call
    poisson_op2_gpu(arg0+n*2,
                arg1+n*1,
                arg2+n*60,
                arg3+n*60,
                arg4+n*60,
                arg5+n*60,
                arg6+n*60,
                arg7+n*60,
                ind_arg0+map8idx*18,
                ind_arg0+map9idx*18,
                ind_arg1+map8idx*1,
                ind_arg1+map9idx*1,
                ind_arg2+map8idx*18,
                ind_arg2+map9idx*18,
                ind_arg3+map8idx*10,
                ind_arg3+map9idx*10,
                arg16_l,
                arg17_l,
                arg18+n*100,
                arg19+n*100);
    atomicAdd(&ind_arg4[0+map8idx*100],arg16_l[0]);
    atomicAdd(&ind_arg4[1+map8idx*100],arg16_l[1]);
    atomicAdd(&ind_arg4[2+map8idx*100],arg16_l[2]);
    atomicAdd(&ind_arg4[3+map8idx*100],arg16_l[3]);
    atomicAdd(&ind_arg4[4+map8idx*100],arg16_l[4]);
    atomicAdd(&ind_arg4[5+map8idx*100],arg16_l[5]);
    atomicAdd(&ind_arg4[6+map8idx*100],arg16_l[6]);
    atomicAdd(&ind_arg4[7+map8idx*100],arg16_l[7]);
    atomicAdd(&ind_arg4[8+map8idx*100],arg16_l[8]);
    atomicAdd(&ind_arg4[9+map8idx*100],arg16_l[9]);
    atomicAdd(&ind_arg4[10+map8idx*100],arg16_l[10]);
    atomicAdd(&ind_arg4[11+map8idx*100],arg16_l[11]);
    atomicAdd(&ind_arg4[12+map8idx*100],arg16_l[12]);
    atomicAdd(&ind_arg4[13+map8idx*100],arg16_l[13]);
    atomicAdd(&ind_arg4[14+map8idx*100],arg16_l[14]);
    atomicAdd(&ind_arg4[15+map8idx*100],arg16_l[15]);
    atomicAdd(&ind_arg4[16+map8idx*100],arg16_l[16]);
    atomicAdd(&ind_arg4[17+map8idx*100],arg16_l[17]);
    atomicAdd(&ind_arg4[18+map8idx*100],arg16_l[18]);
    atomicAdd(&ind_arg4[19+map8idx*100],arg16_l[19]);
    atomicAdd(&ind_arg4[20+map8idx*100],arg16_l[20]);
    atomicAdd(&ind_arg4[21+map8idx*100],arg16_l[21]);
    atomicAdd(&ind_arg4[22+map8idx*100],arg16_l[22]);
    atomicAdd(&ind_arg4[23+map8idx*100],arg16_l[23]);
    atomicAdd(&ind_arg4[24+map8idx*100],arg16_l[24]);
    atomicAdd(&ind_arg4[25+map8idx*100],arg16_l[25]);
    atomicAdd(&ind_arg4[26+map8idx*100],arg16_l[26]);
    atomicAdd(&ind_arg4[27+map8idx*100],arg16_l[27]);
    atomicAdd(&ind_arg4[28+map8idx*100],arg16_l[28]);
    atomicAdd(&ind_arg4[29+map8idx*100],arg16_l[29]);
    atomicAdd(&ind_arg4[30+map8idx*100],arg16_l[30]);
    atomicAdd(&ind_arg4[31+map8idx*100],arg16_l[31]);
    atomicAdd(&ind_arg4[32+map8idx*100],arg16_l[32]);
    atomicAdd(&ind_arg4[33+map8idx*100],arg16_l[33]);
    atomicAdd(&ind_arg4[34+map8idx*100],arg16_l[34]);
    atomicAdd(&ind_arg4[35+map8idx*100],arg16_l[35]);
    atomicAdd(&ind_arg4[36+map8idx*100],arg16_l[36]);
    atomicAdd(&ind_arg4[37+map8idx*100],arg16_l[37]);
    atomicAdd(&ind_arg4[38+map8idx*100],arg16_l[38]);
    atomicAdd(&ind_arg4[39+map8idx*100],arg16_l[39]);
    atomicAdd(&ind_arg4[40+map8idx*100],arg16_l[40]);
    atomicAdd(&ind_arg4[41+map8idx*100],arg16_l[41]);
    atomicAdd(&ind_arg4[42+map8idx*100],arg16_l[42]);
    atomicAdd(&ind_arg4[43+map8idx*100],arg16_l[43]);
    atomicAdd(&ind_arg4[44+map8idx*100],arg16_l[44]);
    atomicAdd(&ind_arg4[45+map8idx*100],arg16_l[45]);
    atomicAdd(&ind_arg4[46+map8idx*100],arg16_l[46]);
    atomicAdd(&ind_arg4[47+map8idx*100],arg16_l[47]);
    atomicAdd(&ind_arg4[48+map8idx*100],arg16_l[48]);
    atomicAdd(&ind_arg4[49+map8idx*100],arg16_l[49]);
    atomicAdd(&ind_arg4[50+map8idx*100],arg16_l[50]);
    atomicAdd(&ind_arg4[51+map8idx*100],arg16_l[51]);
    atomicAdd(&ind_arg4[52+map8idx*100],arg16_l[52]);
    atomicAdd(&ind_arg4[53+map8idx*100],arg16_l[53]);
    atomicAdd(&ind_arg4[54+map8idx*100],arg16_l[54]);
    atomicAdd(&ind_arg4[55+map8idx*100],arg16_l[55]);
    atomicAdd(&ind_arg4[56+map8idx*100],arg16_l[56]);
    atomicAdd(&ind_arg4[57+map8idx*100],arg16_l[57]);
    atomicAdd(&ind_arg4[58+map8idx*100],arg16_l[58]);
    atomicAdd(&ind_arg4[59+map8idx*100],arg16_l[59]);
    atomicAdd(&ind_arg4[60+map8idx*100],arg16_l[60]);
    atomicAdd(&ind_arg4[61+map8idx*100],arg16_l[61]);
    atomicAdd(&ind_arg4[62+map8idx*100],arg16_l[62]);
    atomicAdd(&ind_arg4[63+map8idx*100],arg16_l[63]);
    atomicAdd(&ind_arg4[64+map8idx*100],arg16_l[64]);
    atomicAdd(&ind_arg4[65+map8idx*100],arg16_l[65]);
    atomicAdd(&ind_arg4[66+map8idx*100],arg16_l[66]);
    atomicAdd(&ind_arg4[67+map8idx*100],arg16_l[67]);
    atomicAdd(&ind_arg4[68+map8idx*100],arg16_l[68]);
    atomicAdd(&ind_arg4[69+map8idx*100],arg16_l[69]);
    atomicAdd(&ind_arg4[70+map8idx*100],arg16_l[70]);
    atomicAdd(&ind_arg4[71+map8idx*100],arg16_l[71]);
    atomicAdd(&ind_arg4[72+map8idx*100],arg16_l[72]);
    atomicAdd(&ind_arg4[73+map8idx*100],arg16_l[73]);
    atomicAdd(&ind_arg4[74+map8idx*100],arg16_l[74]);
    atomicAdd(&ind_arg4[75+map8idx*100],arg16_l[75]);
    atomicAdd(&ind_arg4[76+map8idx*100],arg16_l[76]);
    atomicAdd(&ind_arg4[77+map8idx*100],arg16_l[77]);
    atomicAdd(&ind_arg4[78+map8idx*100],arg16_l[78]);
    atomicAdd(&ind_arg4[79+map8idx*100],arg16_l[79]);
    atomicAdd(&ind_arg4[80+map8idx*100],arg16_l[80]);
    atomicAdd(&ind_arg4[81+map8idx*100],arg16_l[81]);
    atomicAdd(&ind_arg4[82+map8idx*100],arg16_l[82]);
    atomicAdd(&ind_arg4[83+map8idx*100],arg16_l[83]);
    atomicAdd(&ind_arg4[84+map8idx*100],arg16_l[84]);
    atomicAdd(&ind_arg4[85+map8idx*100],arg16_l[85]);
    atomicAdd(&ind_arg4[86+map8idx*100],arg16_l[86]);
    atomicAdd(&ind_arg4[87+map8idx*100],arg16_l[87]);
    atomicAdd(&ind_arg4[88+map8idx*100],arg16_l[88]);
    atomicAdd(&ind_arg4[89+map8idx*100],arg16_l[89]);
    atomicAdd(&ind_arg4[90+map8idx*100],arg16_l[90]);
    atomicAdd(&ind_arg4[91+map8idx*100],arg16_l[91]);
    atomicAdd(&ind_arg4[92+map8idx*100],arg16_l[92]);
    atomicAdd(&ind_arg4[93+map8idx*100],arg16_l[93]);
    atomicAdd(&ind_arg4[94+map8idx*100],arg16_l[94]);
    atomicAdd(&ind_arg4[95+map8idx*100],arg16_l[95]);
    atomicAdd(&ind_arg4[96+map8idx*100],arg16_l[96]);
    atomicAdd(&ind_arg4[97+map8idx*100],arg16_l[97]);
    atomicAdd(&ind_arg4[98+map8idx*100],arg16_l[98]);
    atomicAdd(&ind_arg4[99+map8idx*100],arg16_l[99]);
    atomicAdd(&ind_arg4[0+map9idx*100],arg17_l[0]);
    atomicAdd(&ind_arg4[1+map9idx*100],arg17_l[1]);
    atomicAdd(&ind_arg4[2+map9idx*100],arg17_l[2]);
    atomicAdd(&ind_arg4[3+map9idx*100],arg17_l[3]);
    atomicAdd(&ind_arg4[4+map9idx*100],arg17_l[4]);
    atomicAdd(&ind_arg4[5+map9idx*100],arg17_l[5]);
    atomicAdd(&ind_arg4[6+map9idx*100],arg17_l[6]);
    atomicAdd(&ind_arg4[7+map9idx*100],arg17_l[7]);
    atomicAdd(&ind_arg4[8+map9idx*100],arg17_l[8]);
    atomicAdd(&ind_arg4[9+map9idx*100],arg17_l[9]);
    atomicAdd(&ind_arg4[10+map9idx*100],arg17_l[10]);
    atomicAdd(&ind_arg4[11+map9idx*100],arg17_l[11]);
    atomicAdd(&ind_arg4[12+map9idx*100],arg17_l[12]);
    atomicAdd(&ind_arg4[13+map9idx*100],arg17_l[13]);
    atomicAdd(&ind_arg4[14+map9idx*100],arg17_l[14]);
    atomicAdd(&ind_arg4[15+map9idx*100],arg17_l[15]);
    atomicAdd(&ind_arg4[16+map9idx*100],arg17_l[16]);
    atomicAdd(&ind_arg4[17+map9idx*100],arg17_l[17]);
    atomicAdd(&ind_arg4[18+map9idx*100],arg17_l[18]);
    atomicAdd(&ind_arg4[19+map9idx*100],arg17_l[19]);
    atomicAdd(&ind_arg4[20+map9idx*100],arg17_l[20]);
    atomicAdd(&ind_arg4[21+map9idx*100],arg17_l[21]);
    atomicAdd(&ind_arg4[22+map9idx*100],arg17_l[22]);
    atomicAdd(&ind_arg4[23+map9idx*100],arg17_l[23]);
    atomicAdd(&ind_arg4[24+map9idx*100],arg17_l[24]);
    atomicAdd(&ind_arg4[25+map9idx*100],arg17_l[25]);
    atomicAdd(&ind_arg4[26+map9idx*100],arg17_l[26]);
    atomicAdd(&ind_arg4[27+map9idx*100],arg17_l[27]);
    atomicAdd(&ind_arg4[28+map9idx*100],arg17_l[28]);
    atomicAdd(&ind_arg4[29+map9idx*100],arg17_l[29]);
    atomicAdd(&ind_arg4[30+map9idx*100],arg17_l[30]);
    atomicAdd(&ind_arg4[31+map9idx*100],arg17_l[31]);
    atomicAdd(&ind_arg4[32+map9idx*100],arg17_l[32]);
    atomicAdd(&ind_arg4[33+map9idx*100],arg17_l[33]);
    atomicAdd(&ind_arg4[34+map9idx*100],arg17_l[34]);
    atomicAdd(&ind_arg4[35+map9idx*100],arg17_l[35]);
    atomicAdd(&ind_arg4[36+map9idx*100],arg17_l[36]);
    atomicAdd(&ind_arg4[37+map9idx*100],arg17_l[37]);
    atomicAdd(&ind_arg4[38+map9idx*100],arg17_l[38]);
    atomicAdd(&ind_arg4[39+map9idx*100],arg17_l[39]);
    atomicAdd(&ind_arg4[40+map9idx*100],arg17_l[40]);
    atomicAdd(&ind_arg4[41+map9idx*100],arg17_l[41]);
    atomicAdd(&ind_arg4[42+map9idx*100],arg17_l[42]);
    atomicAdd(&ind_arg4[43+map9idx*100],arg17_l[43]);
    atomicAdd(&ind_arg4[44+map9idx*100],arg17_l[44]);
    atomicAdd(&ind_arg4[45+map9idx*100],arg17_l[45]);
    atomicAdd(&ind_arg4[46+map9idx*100],arg17_l[46]);
    atomicAdd(&ind_arg4[47+map9idx*100],arg17_l[47]);
    atomicAdd(&ind_arg4[48+map9idx*100],arg17_l[48]);
    atomicAdd(&ind_arg4[49+map9idx*100],arg17_l[49]);
    atomicAdd(&ind_arg4[50+map9idx*100],arg17_l[50]);
    atomicAdd(&ind_arg4[51+map9idx*100],arg17_l[51]);
    atomicAdd(&ind_arg4[52+map9idx*100],arg17_l[52]);
    atomicAdd(&ind_arg4[53+map9idx*100],arg17_l[53]);
    atomicAdd(&ind_arg4[54+map9idx*100],arg17_l[54]);
    atomicAdd(&ind_arg4[55+map9idx*100],arg17_l[55]);
    atomicAdd(&ind_arg4[56+map9idx*100],arg17_l[56]);
    atomicAdd(&ind_arg4[57+map9idx*100],arg17_l[57]);
    atomicAdd(&ind_arg4[58+map9idx*100],arg17_l[58]);
    atomicAdd(&ind_arg4[59+map9idx*100],arg17_l[59]);
    atomicAdd(&ind_arg4[60+map9idx*100],arg17_l[60]);
    atomicAdd(&ind_arg4[61+map9idx*100],arg17_l[61]);
    atomicAdd(&ind_arg4[62+map9idx*100],arg17_l[62]);
    atomicAdd(&ind_arg4[63+map9idx*100],arg17_l[63]);
    atomicAdd(&ind_arg4[64+map9idx*100],arg17_l[64]);
    atomicAdd(&ind_arg4[65+map9idx*100],arg17_l[65]);
    atomicAdd(&ind_arg4[66+map9idx*100],arg17_l[66]);
    atomicAdd(&ind_arg4[67+map9idx*100],arg17_l[67]);
    atomicAdd(&ind_arg4[68+map9idx*100],arg17_l[68]);
    atomicAdd(&ind_arg4[69+map9idx*100],arg17_l[69]);
    atomicAdd(&ind_arg4[70+map9idx*100],arg17_l[70]);
    atomicAdd(&ind_arg4[71+map9idx*100],arg17_l[71]);
    atomicAdd(&ind_arg4[72+map9idx*100],arg17_l[72]);
    atomicAdd(&ind_arg4[73+map9idx*100],arg17_l[73]);
    atomicAdd(&ind_arg4[74+map9idx*100],arg17_l[74]);
    atomicAdd(&ind_arg4[75+map9idx*100],arg17_l[75]);
    atomicAdd(&ind_arg4[76+map9idx*100],arg17_l[76]);
    atomicAdd(&ind_arg4[77+map9idx*100],arg17_l[77]);
    atomicAdd(&ind_arg4[78+map9idx*100],arg17_l[78]);
    atomicAdd(&ind_arg4[79+map9idx*100],arg17_l[79]);
    atomicAdd(&ind_arg4[80+map9idx*100],arg17_l[80]);
    atomicAdd(&ind_arg4[81+map9idx*100],arg17_l[81]);
    atomicAdd(&ind_arg4[82+map9idx*100],arg17_l[82]);
    atomicAdd(&ind_arg4[83+map9idx*100],arg17_l[83]);
    atomicAdd(&ind_arg4[84+map9idx*100],arg17_l[84]);
    atomicAdd(&ind_arg4[85+map9idx*100],arg17_l[85]);
    atomicAdd(&ind_arg4[86+map9idx*100],arg17_l[86]);
    atomicAdd(&ind_arg4[87+map9idx*100],arg17_l[87]);
    atomicAdd(&ind_arg4[88+map9idx*100],arg17_l[88]);
    atomicAdd(&ind_arg4[89+map9idx*100],arg17_l[89]);
    atomicAdd(&ind_arg4[90+map9idx*100],arg17_l[90]);
    atomicAdd(&ind_arg4[91+map9idx*100],arg17_l[91]);
    atomicAdd(&ind_arg4[92+map9idx*100],arg17_l[92]);
    atomicAdd(&ind_arg4[93+map9idx*100],arg17_l[93]);
    atomicAdd(&ind_arg4[94+map9idx*100],arg17_l[94]);
    atomicAdd(&ind_arg4[95+map9idx*100],arg17_l[95]);
    atomicAdd(&ind_arg4[96+map9idx*100],arg17_l[96]);
    atomicAdd(&ind_arg4[97+map9idx*100],arg17_l[97]);
    atomicAdd(&ind_arg4[98+map9idx*100],arg17_l[98]);
    atomicAdd(&ind_arg4[99+map9idx*100],arg17_l[99]);
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


  int    ninds   = 5;
  int    inds[20] = {-1,-1,-1,-1,-1,-1,-1,-1,0,0,1,1,2,2,3,3,4,4,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_op2\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_23
      int nthread = OP_BLOCK_SIZE_23;
    #else
      int nthread = OP_block_size;
    #endif

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_grouped(nargs, args, 2);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        int nblocks = (end-start-1)/nthread+1;
        op_cuda_poisson_op2<<<nblocks,nthread>>>(
        (double *)arg8.data_d,
        (double *)arg10.data_d,
        (double *)arg12.data_d,
        (double *)arg14.data_d,
        (double *)arg16.data_d,
        arg8.map_data_d,
        (int*)arg0.data_d,
        (bool*)arg1.data_d,
        (double*)arg2.data_d,
        (double*)arg3.data_d,
        (double*)arg4.data_d,
        (double*)arg5.data_d,
        (double*)arg6.data_d,
        (double*)arg7.data_d,
        (double*)arg18.data_d,
        (double*)arg19.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[23].time     += wall_t2 - wall_t1;
}
