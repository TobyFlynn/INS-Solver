//
// auto-generated by op2.py
//

//user function
__device__ void gauss_gfi_faces_gpu( const int *edgeNum, const bool *rev,
                            double **gf0, double **gf1, double **gf2) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int indL, indR;
      if(!reverse) {
        indL = m * 10 + n;
        indR = m * 10 + n;
      } else {
        indL = m * 10 + n;
        indR = (6 - 1 - m) * 10 + n;
      }

      if(edgeL == 0) {
        if(edgeR == 0) {
          gf0[0][indL] += gFInterp0_g_cuda[indR];
          gf0[1][indR] += gFInterp0_g_cuda[indL];
        } else if(edgeR == 1) {
          gf0[0][indL] += gFInterp1_g_cuda[indR];
          gf1[1][indR] += gFInterp0_g_cuda[indL];
        } else {
          gf0[0][indL] += gFInterp2_g_cuda[indR];
          gf2[1][indR] += gFInterp0_g_cuda[indL];
        }
      } else if(edgeL == 1) {
        if(edgeR == 0) {
          gf1[0][indL] += gFInterp0_g_cuda[indR];
          gf0[1][indR] += gFInterp1_g_cuda[indL];
        } else if(edgeR == 1) {
          gf1[0][indL] += gFInterp1_g_cuda[indR];
          gf1[1][indR] += gFInterp1_g_cuda[indL];
        } else {
          gf1[0][indL] += gFInterp2_g_cuda[indR];
          gf2[1][indR] += gFInterp1_g_cuda[indL];
        }
      } else {
        if(edgeR == 0) {
          gf2[0][indL] += gFInterp0_g_cuda[indR];
          gf0[1][indR] += gFInterp2_g_cuda[indL];
        } else if(edgeR == 1) {
          gf2[0][indL] += gFInterp1_g_cuda[indR];
          gf1[1][indR] += gFInterp2_g_cuda[indL];
        } else {
          gf2[0][indL] += gFInterp2_g_cuda[indR];
          gf2[1][indR] += gFInterp2_g_cuda[indL];
        }
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_gauss_gfi_faces(
  double *__restrict ind_arg0,
  double *__restrict ind_arg1,
  double *__restrict ind_arg2,
  const int *__restrict opDat2Map,
  const int *__restrict arg0,
  const bool *__restrict arg1,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg2_l[60];
    for ( int d=0; d<60; d++ ){
      arg2_l[d] = ZERO_double;
    }
    double arg3_l[60];
    for ( int d=0; d<60; d++ ){
      arg3_l[d] = ZERO_double;
    }
    double arg4_l[60];
    for ( int d=0; d<60; d++ ){
      arg4_l[d] = ZERO_double;
    }
    double arg5_l[60];
    for ( int d=0; d<60; d++ ){
      arg5_l[d] = ZERO_double;
    }
    double arg6_l[60];
    for ( int d=0; d<60; d++ ){
      arg6_l[d] = ZERO_double;
    }
    double arg7_l[60];
    for ( int d=0; d<60; d++ ){
      arg7_l[d] = ZERO_double;
    }
    int map2idx;
    int map3idx;
    map2idx = opDat2Map[n + set_size * 0];
    map3idx = opDat2Map[n + set_size * 1];
    double* arg2_vec[] = {
      arg2_l,
      arg3_l};
    double* arg4_vec[] = {
      arg4_l,
      arg5_l};
    double* arg6_vec[] = {
      arg6_l,
      arg7_l};

    //user-supplied kernel call
    gauss_gfi_faces_gpu(arg0+n*2,
                    arg1+n*1,
                    arg2_vec,
                    arg4_vec,
                    arg6_vec);
    atomicAdd(&ind_arg0[0+map2idx*60],arg2_l[0]);
    atomicAdd(&ind_arg0[1+map2idx*60],arg2_l[1]);
    atomicAdd(&ind_arg0[2+map2idx*60],arg2_l[2]);
    atomicAdd(&ind_arg0[3+map2idx*60],arg2_l[3]);
    atomicAdd(&ind_arg0[4+map2idx*60],arg2_l[4]);
    atomicAdd(&ind_arg0[5+map2idx*60],arg2_l[5]);
    atomicAdd(&ind_arg0[6+map2idx*60],arg2_l[6]);
    atomicAdd(&ind_arg0[7+map2idx*60],arg2_l[7]);
    atomicAdd(&ind_arg0[8+map2idx*60],arg2_l[8]);
    atomicAdd(&ind_arg0[9+map2idx*60],arg2_l[9]);
    atomicAdd(&ind_arg0[10+map2idx*60],arg2_l[10]);
    atomicAdd(&ind_arg0[11+map2idx*60],arg2_l[11]);
    atomicAdd(&ind_arg0[12+map2idx*60],arg2_l[12]);
    atomicAdd(&ind_arg0[13+map2idx*60],arg2_l[13]);
    atomicAdd(&ind_arg0[14+map2idx*60],arg2_l[14]);
    atomicAdd(&ind_arg0[15+map2idx*60],arg2_l[15]);
    atomicAdd(&ind_arg0[16+map2idx*60],arg2_l[16]);
    atomicAdd(&ind_arg0[17+map2idx*60],arg2_l[17]);
    atomicAdd(&ind_arg0[18+map2idx*60],arg2_l[18]);
    atomicAdd(&ind_arg0[19+map2idx*60],arg2_l[19]);
    atomicAdd(&ind_arg0[20+map2idx*60],arg2_l[20]);
    atomicAdd(&ind_arg0[21+map2idx*60],arg2_l[21]);
    atomicAdd(&ind_arg0[22+map2idx*60],arg2_l[22]);
    atomicAdd(&ind_arg0[23+map2idx*60],arg2_l[23]);
    atomicAdd(&ind_arg0[24+map2idx*60],arg2_l[24]);
    atomicAdd(&ind_arg0[25+map2idx*60],arg2_l[25]);
    atomicAdd(&ind_arg0[26+map2idx*60],arg2_l[26]);
    atomicAdd(&ind_arg0[27+map2idx*60],arg2_l[27]);
    atomicAdd(&ind_arg0[28+map2idx*60],arg2_l[28]);
    atomicAdd(&ind_arg0[29+map2idx*60],arg2_l[29]);
    atomicAdd(&ind_arg0[30+map2idx*60],arg2_l[30]);
    atomicAdd(&ind_arg0[31+map2idx*60],arg2_l[31]);
    atomicAdd(&ind_arg0[32+map2idx*60],arg2_l[32]);
    atomicAdd(&ind_arg0[33+map2idx*60],arg2_l[33]);
    atomicAdd(&ind_arg0[34+map2idx*60],arg2_l[34]);
    atomicAdd(&ind_arg0[35+map2idx*60],arg2_l[35]);
    atomicAdd(&ind_arg0[36+map2idx*60],arg2_l[36]);
    atomicAdd(&ind_arg0[37+map2idx*60],arg2_l[37]);
    atomicAdd(&ind_arg0[38+map2idx*60],arg2_l[38]);
    atomicAdd(&ind_arg0[39+map2idx*60],arg2_l[39]);
    atomicAdd(&ind_arg0[40+map2idx*60],arg2_l[40]);
    atomicAdd(&ind_arg0[41+map2idx*60],arg2_l[41]);
    atomicAdd(&ind_arg0[42+map2idx*60],arg2_l[42]);
    atomicAdd(&ind_arg0[43+map2idx*60],arg2_l[43]);
    atomicAdd(&ind_arg0[44+map2idx*60],arg2_l[44]);
    atomicAdd(&ind_arg0[45+map2idx*60],arg2_l[45]);
    atomicAdd(&ind_arg0[46+map2idx*60],arg2_l[46]);
    atomicAdd(&ind_arg0[47+map2idx*60],arg2_l[47]);
    atomicAdd(&ind_arg0[48+map2idx*60],arg2_l[48]);
    atomicAdd(&ind_arg0[49+map2idx*60],arg2_l[49]);
    atomicAdd(&ind_arg0[50+map2idx*60],arg2_l[50]);
    atomicAdd(&ind_arg0[51+map2idx*60],arg2_l[51]);
    atomicAdd(&ind_arg0[52+map2idx*60],arg2_l[52]);
    atomicAdd(&ind_arg0[53+map2idx*60],arg2_l[53]);
    atomicAdd(&ind_arg0[54+map2idx*60],arg2_l[54]);
    atomicAdd(&ind_arg0[55+map2idx*60],arg2_l[55]);
    atomicAdd(&ind_arg0[56+map2idx*60],arg2_l[56]);
    atomicAdd(&ind_arg0[57+map2idx*60],arg2_l[57]);
    atomicAdd(&ind_arg0[58+map2idx*60],arg2_l[58]);
    atomicAdd(&ind_arg0[59+map2idx*60],arg2_l[59]);
    atomicAdd(&ind_arg0[0+map3idx*60],arg3_l[0]);
    atomicAdd(&ind_arg0[1+map3idx*60],arg3_l[1]);
    atomicAdd(&ind_arg0[2+map3idx*60],arg3_l[2]);
    atomicAdd(&ind_arg0[3+map3idx*60],arg3_l[3]);
    atomicAdd(&ind_arg0[4+map3idx*60],arg3_l[4]);
    atomicAdd(&ind_arg0[5+map3idx*60],arg3_l[5]);
    atomicAdd(&ind_arg0[6+map3idx*60],arg3_l[6]);
    atomicAdd(&ind_arg0[7+map3idx*60],arg3_l[7]);
    atomicAdd(&ind_arg0[8+map3idx*60],arg3_l[8]);
    atomicAdd(&ind_arg0[9+map3idx*60],arg3_l[9]);
    atomicAdd(&ind_arg0[10+map3idx*60],arg3_l[10]);
    atomicAdd(&ind_arg0[11+map3idx*60],arg3_l[11]);
    atomicAdd(&ind_arg0[12+map3idx*60],arg3_l[12]);
    atomicAdd(&ind_arg0[13+map3idx*60],arg3_l[13]);
    atomicAdd(&ind_arg0[14+map3idx*60],arg3_l[14]);
    atomicAdd(&ind_arg0[15+map3idx*60],arg3_l[15]);
    atomicAdd(&ind_arg0[16+map3idx*60],arg3_l[16]);
    atomicAdd(&ind_arg0[17+map3idx*60],arg3_l[17]);
    atomicAdd(&ind_arg0[18+map3idx*60],arg3_l[18]);
    atomicAdd(&ind_arg0[19+map3idx*60],arg3_l[19]);
    atomicAdd(&ind_arg0[20+map3idx*60],arg3_l[20]);
    atomicAdd(&ind_arg0[21+map3idx*60],arg3_l[21]);
    atomicAdd(&ind_arg0[22+map3idx*60],arg3_l[22]);
    atomicAdd(&ind_arg0[23+map3idx*60],arg3_l[23]);
    atomicAdd(&ind_arg0[24+map3idx*60],arg3_l[24]);
    atomicAdd(&ind_arg0[25+map3idx*60],arg3_l[25]);
    atomicAdd(&ind_arg0[26+map3idx*60],arg3_l[26]);
    atomicAdd(&ind_arg0[27+map3idx*60],arg3_l[27]);
    atomicAdd(&ind_arg0[28+map3idx*60],arg3_l[28]);
    atomicAdd(&ind_arg0[29+map3idx*60],arg3_l[29]);
    atomicAdd(&ind_arg0[30+map3idx*60],arg3_l[30]);
    atomicAdd(&ind_arg0[31+map3idx*60],arg3_l[31]);
    atomicAdd(&ind_arg0[32+map3idx*60],arg3_l[32]);
    atomicAdd(&ind_arg0[33+map3idx*60],arg3_l[33]);
    atomicAdd(&ind_arg0[34+map3idx*60],arg3_l[34]);
    atomicAdd(&ind_arg0[35+map3idx*60],arg3_l[35]);
    atomicAdd(&ind_arg0[36+map3idx*60],arg3_l[36]);
    atomicAdd(&ind_arg0[37+map3idx*60],arg3_l[37]);
    atomicAdd(&ind_arg0[38+map3idx*60],arg3_l[38]);
    atomicAdd(&ind_arg0[39+map3idx*60],arg3_l[39]);
    atomicAdd(&ind_arg0[40+map3idx*60],arg3_l[40]);
    atomicAdd(&ind_arg0[41+map3idx*60],arg3_l[41]);
    atomicAdd(&ind_arg0[42+map3idx*60],arg3_l[42]);
    atomicAdd(&ind_arg0[43+map3idx*60],arg3_l[43]);
    atomicAdd(&ind_arg0[44+map3idx*60],arg3_l[44]);
    atomicAdd(&ind_arg0[45+map3idx*60],arg3_l[45]);
    atomicAdd(&ind_arg0[46+map3idx*60],arg3_l[46]);
    atomicAdd(&ind_arg0[47+map3idx*60],arg3_l[47]);
    atomicAdd(&ind_arg0[48+map3idx*60],arg3_l[48]);
    atomicAdd(&ind_arg0[49+map3idx*60],arg3_l[49]);
    atomicAdd(&ind_arg0[50+map3idx*60],arg3_l[50]);
    atomicAdd(&ind_arg0[51+map3idx*60],arg3_l[51]);
    atomicAdd(&ind_arg0[52+map3idx*60],arg3_l[52]);
    atomicAdd(&ind_arg0[53+map3idx*60],arg3_l[53]);
    atomicAdd(&ind_arg0[54+map3idx*60],arg3_l[54]);
    atomicAdd(&ind_arg0[55+map3idx*60],arg3_l[55]);
    atomicAdd(&ind_arg0[56+map3idx*60],arg3_l[56]);
    atomicAdd(&ind_arg0[57+map3idx*60],arg3_l[57]);
    atomicAdd(&ind_arg0[58+map3idx*60],arg3_l[58]);
    atomicAdd(&ind_arg0[59+map3idx*60],arg3_l[59]);
    atomicAdd(&ind_arg1[0+map2idx*60],arg4_l[0]);
    atomicAdd(&ind_arg1[1+map2idx*60],arg4_l[1]);
    atomicAdd(&ind_arg1[2+map2idx*60],arg4_l[2]);
    atomicAdd(&ind_arg1[3+map2idx*60],arg4_l[3]);
    atomicAdd(&ind_arg1[4+map2idx*60],arg4_l[4]);
    atomicAdd(&ind_arg1[5+map2idx*60],arg4_l[5]);
    atomicAdd(&ind_arg1[6+map2idx*60],arg4_l[6]);
    atomicAdd(&ind_arg1[7+map2idx*60],arg4_l[7]);
    atomicAdd(&ind_arg1[8+map2idx*60],arg4_l[8]);
    atomicAdd(&ind_arg1[9+map2idx*60],arg4_l[9]);
    atomicAdd(&ind_arg1[10+map2idx*60],arg4_l[10]);
    atomicAdd(&ind_arg1[11+map2idx*60],arg4_l[11]);
    atomicAdd(&ind_arg1[12+map2idx*60],arg4_l[12]);
    atomicAdd(&ind_arg1[13+map2idx*60],arg4_l[13]);
    atomicAdd(&ind_arg1[14+map2idx*60],arg4_l[14]);
    atomicAdd(&ind_arg1[15+map2idx*60],arg4_l[15]);
    atomicAdd(&ind_arg1[16+map2idx*60],arg4_l[16]);
    atomicAdd(&ind_arg1[17+map2idx*60],arg4_l[17]);
    atomicAdd(&ind_arg1[18+map2idx*60],arg4_l[18]);
    atomicAdd(&ind_arg1[19+map2idx*60],arg4_l[19]);
    atomicAdd(&ind_arg1[20+map2idx*60],arg4_l[20]);
    atomicAdd(&ind_arg1[21+map2idx*60],arg4_l[21]);
    atomicAdd(&ind_arg1[22+map2idx*60],arg4_l[22]);
    atomicAdd(&ind_arg1[23+map2idx*60],arg4_l[23]);
    atomicAdd(&ind_arg1[24+map2idx*60],arg4_l[24]);
    atomicAdd(&ind_arg1[25+map2idx*60],arg4_l[25]);
    atomicAdd(&ind_arg1[26+map2idx*60],arg4_l[26]);
    atomicAdd(&ind_arg1[27+map2idx*60],arg4_l[27]);
    atomicAdd(&ind_arg1[28+map2idx*60],arg4_l[28]);
    atomicAdd(&ind_arg1[29+map2idx*60],arg4_l[29]);
    atomicAdd(&ind_arg1[30+map2idx*60],arg4_l[30]);
    atomicAdd(&ind_arg1[31+map2idx*60],arg4_l[31]);
    atomicAdd(&ind_arg1[32+map2idx*60],arg4_l[32]);
    atomicAdd(&ind_arg1[33+map2idx*60],arg4_l[33]);
    atomicAdd(&ind_arg1[34+map2idx*60],arg4_l[34]);
    atomicAdd(&ind_arg1[35+map2idx*60],arg4_l[35]);
    atomicAdd(&ind_arg1[36+map2idx*60],arg4_l[36]);
    atomicAdd(&ind_arg1[37+map2idx*60],arg4_l[37]);
    atomicAdd(&ind_arg1[38+map2idx*60],arg4_l[38]);
    atomicAdd(&ind_arg1[39+map2idx*60],arg4_l[39]);
    atomicAdd(&ind_arg1[40+map2idx*60],arg4_l[40]);
    atomicAdd(&ind_arg1[41+map2idx*60],arg4_l[41]);
    atomicAdd(&ind_arg1[42+map2idx*60],arg4_l[42]);
    atomicAdd(&ind_arg1[43+map2idx*60],arg4_l[43]);
    atomicAdd(&ind_arg1[44+map2idx*60],arg4_l[44]);
    atomicAdd(&ind_arg1[45+map2idx*60],arg4_l[45]);
    atomicAdd(&ind_arg1[46+map2idx*60],arg4_l[46]);
    atomicAdd(&ind_arg1[47+map2idx*60],arg4_l[47]);
    atomicAdd(&ind_arg1[48+map2idx*60],arg4_l[48]);
    atomicAdd(&ind_arg1[49+map2idx*60],arg4_l[49]);
    atomicAdd(&ind_arg1[50+map2idx*60],arg4_l[50]);
    atomicAdd(&ind_arg1[51+map2idx*60],arg4_l[51]);
    atomicAdd(&ind_arg1[52+map2idx*60],arg4_l[52]);
    atomicAdd(&ind_arg1[53+map2idx*60],arg4_l[53]);
    atomicAdd(&ind_arg1[54+map2idx*60],arg4_l[54]);
    atomicAdd(&ind_arg1[55+map2idx*60],arg4_l[55]);
    atomicAdd(&ind_arg1[56+map2idx*60],arg4_l[56]);
    atomicAdd(&ind_arg1[57+map2idx*60],arg4_l[57]);
    atomicAdd(&ind_arg1[58+map2idx*60],arg4_l[58]);
    atomicAdd(&ind_arg1[59+map2idx*60],arg4_l[59]);
    atomicAdd(&ind_arg1[0+map3idx*60],arg5_l[0]);
    atomicAdd(&ind_arg1[1+map3idx*60],arg5_l[1]);
    atomicAdd(&ind_arg1[2+map3idx*60],arg5_l[2]);
    atomicAdd(&ind_arg1[3+map3idx*60],arg5_l[3]);
    atomicAdd(&ind_arg1[4+map3idx*60],arg5_l[4]);
    atomicAdd(&ind_arg1[5+map3idx*60],arg5_l[5]);
    atomicAdd(&ind_arg1[6+map3idx*60],arg5_l[6]);
    atomicAdd(&ind_arg1[7+map3idx*60],arg5_l[7]);
    atomicAdd(&ind_arg1[8+map3idx*60],arg5_l[8]);
    atomicAdd(&ind_arg1[9+map3idx*60],arg5_l[9]);
    atomicAdd(&ind_arg1[10+map3idx*60],arg5_l[10]);
    atomicAdd(&ind_arg1[11+map3idx*60],arg5_l[11]);
    atomicAdd(&ind_arg1[12+map3idx*60],arg5_l[12]);
    atomicAdd(&ind_arg1[13+map3idx*60],arg5_l[13]);
    atomicAdd(&ind_arg1[14+map3idx*60],arg5_l[14]);
    atomicAdd(&ind_arg1[15+map3idx*60],arg5_l[15]);
    atomicAdd(&ind_arg1[16+map3idx*60],arg5_l[16]);
    atomicAdd(&ind_arg1[17+map3idx*60],arg5_l[17]);
    atomicAdd(&ind_arg1[18+map3idx*60],arg5_l[18]);
    atomicAdd(&ind_arg1[19+map3idx*60],arg5_l[19]);
    atomicAdd(&ind_arg1[20+map3idx*60],arg5_l[20]);
    atomicAdd(&ind_arg1[21+map3idx*60],arg5_l[21]);
    atomicAdd(&ind_arg1[22+map3idx*60],arg5_l[22]);
    atomicAdd(&ind_arg1[23+map3idx*60],arg5_l[23]);
    atomicAdd(&ind_arg1[24+map3idx*60],arg5_l[24]);
    atomicAdd(&ind_arg1[25+map3idx*60],arg5_l[25]);
    atomicAdd(&ind_arg1[26+map3idx*60],arg5_l[26]);
    atomicAdd(&ind_arg1[27+map3idx*60],arg5_l[27]);
    atomicAdd(&ind_arg1[28+map3idx*60],arg5_l[28]);
    atomicAdd(&ind_arg1[29+map3idx*60],arg5_l[29]);
    atomicAdd(&ind_arg1[30+map3idx*60],arg5_l[30]);
    atomicAdd(&ind_arg1[31+map3idx*60],arg5_l[31]);
    atomicAdd(&ind_arg1[32+map3idx*60],arg5_l[32]);
    atomicAdd(&ind_arg1[33+map3idx*60],arg5_l[33]);
    atomicAdd(&ind_arg1[34+map3idx*60],arg5_l[34]);
    atomicAdd(&ind_arg1[35+map3idx*60],arg5_l[35]);
    atomicAdd(&ind_arg1[36+map3idx*60],arg5_l[36]);
    atomicAdd(&ind_arg1[37+map3idx*60],arg5_l[37]);
    atomicAdd(&ind_arg1[38+map3idx*60],arg5_l[38]);
    atomicAdd(&ind_arg1[39+map3idx*60],arg5_l[39]);
    atomicAdd(&ind_arg1[40+map3idx*60],arg5_l[40]);
    atomicAdd(&ind_arg1[41+map3idx*60],arg5_l[41]);
    atomicAdd(&ind_arg1[42+map3idx*60],arg5_l[42]);
    atomicAdd(&ind_arg1[43+map3idx*60],arg5_l[43]);
    atomicAdd(&ind_arg1[44+map3idx*60],arg5_l[44]);
    atomicAdd(&ind_arg1[45+map3idx*60],arg5_l[45]);
    atomicAdd(&ind_arg1[46+map3idx*60],arg5_l[46]);
    atomicAdd(&ind_arg1[47+map3idx*60],arg5_l[47]);
    atomicAdd(&ind_arg1[48+map3idx*60],arg5_l[48]);
    atomicAdd(&ind_arg1[49+map3idx*60],arg5_l[49]);
    atomicAdd(&ind_arg1[50+map3idx*60],arg5_l[50]);
    atomicAdd(&ind_arg1[51+map3idx*60],arg5_l[51]);
    atomicAdd(&ind_arg1[52+map3idx*60],arg5_l[52]);
    atomicAdd(&ind_arg1[53+map3idx*60],arg5_l[53]);
    atomicAdd(&ind_arg1[54+map3idx*60],arg5_l[54]);
    atomicAdd(&ind_arg1[55+map3idx*60],arg5_l[55]);
    atomicAdd(&ind_arg1[56+map3idx*60],arg5_l[56]);
    atomicAdd(&ind_arg1[57+map3idx*60],arg5_l[57]);
    atomicAdd(&ind_arg1[58+map3idx*60],arg5_l[58]);
    atomicAdd(&ind_arg1[59+map3idx*60],arg5_l[59]);
    atomicAdd(&ind_arg2[0+map2idx*60],arg6_l[0]);
    atomicAdd(&ind_arg2[1+map2idx*60],arg6_l[1]);
    atomicAdd(&ind_arg2[2+map2idx*60],arg6_l[2]);
    atomicAdd(&ind_arg2[3+map2idx*60],arg6_l[3]);
    atomicAdd(&ind_arg2[4+map2idx*60],arg6_l[4]);
    atomicAdd(&ind_arg2[5+map2idx*60],arg6_l[5]);
    atomicAdd(&ind_arg2[6+map2idx*60],arg6_l[6]);
    atomicAdd(&ind_arg2[7+map2idx*60],arg6_l[7]);
    atomicAdd(&ind_arg2[8+map2idx*60],arg6_l[8]);
    atomicAdd(&ind_arg2[9+map2idx*60],arg6_l[9]);
    atomicAdd(&ind_arg2[10+map2idx*60],arg6_l[10]);
    atomicAdd(&ind_arg2[11+map2idx*60],arg6_l[11]);
    atomicAdd(&ind_arg2[12+map2idx*60],arg6_l[12]);
    atomicAdd(&ind_arg2[13+map2idx*60],arg6_l[13]);
    atomicAdd(&ind_arg2[14+map2idx*60],arg6_l[14]);
    atomicAdd(&ind_arg2[15+map2idx*60],arg6_l[15]);
    atomicAdd(&ind_arg2[16+map2idx*60],arg6_l[16]);
    atomicAdd(&ind_arg2[17+map2idx*60],arg6_l[17]);
    atomicAdd(&ind_arg2[18+map2idx*60],arg6_l[18]);
    atomicAdd(&ind_arg2[19+map2idx*60],arg6_l[19]);
    atomicAdd(&ind_arg2[20+map2idx*60],arg6_l[20]);
    atomicAdd(&ind_arg2[21+map2idx*60],arg6_l[21]);
    atomicAdd(&ind_arg2[22+map2idx*60],arg6_l[22]);
    atomicAdd(&ind_arg2[23+map2idx*60],arg6_l[23]);
    atomicAdd(&ind_arg2[24+map2idx*60],arg6_l[24]);
    atomicAdd(&ind_arg2[25+map2idx*60],arg6_l[25]);
    atomicAdd(&ind_arg2[26+map2idx*60],arg6_l[26]);
    atomicAdd(&ind_arg2[27+map2idx*60],arg6_l[27]);
    atomicAdd(&ind_arg2[28+map2idx*60],arg6_l[28]);
    atomicAdd(&ind_arg2[29+map2idx*60],arg6_l[29]);
    atomicAdd(&ind_arg2[30+map2idx*60],arg6_l[30]);
    atomicAdd(&ind_arg2[31+map2idx*60],arg6_l[31]);
    atomicAdd(&ind_arg2[32+map2idx*60],arg6_l[32]);
    atomicAdd(&ind_arg2[33+map2idx*60],arg6_l[33]);
    atomicAdd(&ind_arg2[34+map2idx*60],arg6_l[34]);
    atomicAdd(&ind_arg2[35+map2idx*60],arg6_l[35]);
    atomicAdd(&ind_arg2[36+map2idx*60],arg6_l[36]);
    atomicAdd(&ind_arg2[37+map2idx*60],arg6_l[37]);
    atomicAdd(&ind_arg2[38+map2idx*60],arg6_l[38]);
    atomicAdd(&ind_arg2[39+map2idx*60],arg6_l[39]);
    atomicAdd(&ind_arg2[40+map2idx*60],arg6_l[40]);
    atomicAdd(&ind_arg2[41+map2idx*60],arg6_l[41]);
    atomicAdd(&ind_arg2[42+map2idx*60],arg6_l[42]);
    atomicAdd(&ind_arg2[43+map2idx*60],arg6_l[43]);
    atomicAdd(&ind_arg2[44+map2idx*60],arg6_l[44]);
    atomicAdd(&ind_arg2[45+map2idx*60],arg6_l[45]);
    atomicAdd(&ind_arg2[46+map2idx*60],arg6_l[46]);
    atomicAdd(&ind_arg2[47+map2idx*60],arg6_l[47]);
    atomicAdd(&ind_arg2[48+map2idx*60],arg6_l[48]);
    atomicAdd(&ind_arg2[49+map2idx*60],arg6_l[49]);
    atomicAdd(&ind_arg2[50+map2idx*60],arg6_l[50]);
    atomicAdd(&ind_arg2[51+map2idx*60],arg6_l[51]);
    atomicAdd(&ind_arg2[52+map2idx*60],arg6_l[52]);
    atomicAdd(&ind_arg2[53+map2idx*60],arg6_l[53]);
    atomicAdd(&ind_arg2[54+map2idx*60],arg6_l[54]);
    atomicAdd(&ind_arg2[55+map2idx*60],arg6_l[55]);
    atomicAdd(&ind_arg2[56+map2idx*60],arg6_l[56]);
    atomicAdd(&ind_arg2[57+map2idx*60],arg6_l[57]);
    atomicAdd(&ind_arg2[58+map2idx*60],arg6_l[58]);
    atomicAdd(&ind_arg2[59+map2idx*60],arg6_l[59]);
    atomicAdd(&ind_arg2[0+map3idx*60],arg7_l[0]);
    atomicAdd(&ind_arg2[1+map3idx*60],arg7_l[1]);
    atomicAdd(&ind_arg2[2+map3idx*60],arg7_l[2]);
    atomicAdd(&ind_arg2[3+map3idx*60],arg7_l[3]);
    atomicAdd(&ind_arg2[4+map3idx*60],arg7_l[4]);
    atomicAdd(&ind_arg2[5+map3idx*60],arg7_l[5]);
    atomicAdd(&ind_arg2[6+map3idx*60],arg7_l[6]);
    atomicAdd(&ind_arg2[7+map3idx*60],arg7_l[7]);
    atomicAdd(&ind_arg2[8+map3idx*60],arg7_l[8]);
    atomicAdd(&ind_arg2[9+map3idx*60],arg7_l[9]);
    atomicAdd(&ind_arg2[10+map3idx*60],arg7_l[10]);
    atomicAdd(&ind_arg2[11+map3idx*60],arg7_l[11]);
    atomicAdd(&ind_arg2[12+map3idx*60],arg7_l[12]);
    atomicAdd(&ind_arg2[13+map3idx*60],arg7_l[13]);
    atomicAdd(&ind_arg2[14+map3idx*60],arg7_l[14]);
    atomicAdd(&ind_arg2[15+map3idx*60],arg7_l[15]);
    atomicAdd(&ind_arg2[16+map3idx*60],arg7_l[16]);
    atomicAdd(&ind_arg2[17+map3idx*60],arg7_l[17]);
    atomicAdd(&ind_arg2[18+map3idx*60],arg7_l[18]);
    atomicAdd(&ind_arg2[19+map3idx*60],arg7_l[19]);
    atomicAdd(&ind_arg2[20+map3idx*60],arg7_l[20]);
    atomicAdd(&ind_arg2[21+map3idx*60],arg7_l[21]);
    atomicAdd(&ind_arg2[22+map3idx*60],arg7_l[22]);
    atomicAdd(&ind_arg2[23+map3idx*60],arg7_l[23]);
    atomicAdd(&ind_arg2[24+map3idx*60],arg7_l[24]);
    atomicAdd(&ind_arg2[25+map3idx*60],arg7_l[25]);
    atomicAdd(&ind_arg2[26+map3idx*60],arg7_l[26]);
    atomicAdd(&ind_arg2[27+map3idx*60],arg7_l[27]);
    atomicAdd(&ind_arg2[28+map3idx*60],arg7_l[28]);
    atomicAdd(&ind_arg2[29+map3idx*60],arg7_l[29]);
    atomicAdd(&ind_arg2[30+map3idx*60],arg7_l[30]);
    atomicAdd(&ind_arg2[31+map3idx*60],arg7_l[31]);
    atomicAdd(&ind_arg2[32+map3idx*60],arg7_l[32]);
    atomicAdd(&ind_arg2[33+map3idx*60],arg7_l[33]);
    atomicAdd(&ind_arg2[34+map3idx*60],arg7_l[34]);
    atomicAdd(&ind_arg2[35+map3idx*60],arg7_l[35]);
    atomicAdd(&ind_arg2[36+map3idx*60],arg7_l[36]);
    atomicAdd(&ind_arg2[37+map3idx*60],arg7_l[37]);
    atomicAdd(&ind_arg2[38+map3idx*60],arg7_l[38]);
    atomicAdd(&ind_arg2[39+map3idx*60],arg7_l[39]);
    atomicAdd(&ind_arg2[40+map3idx*60],arg7_l[40]);
    atomicAdd(&ind_arg2[41+map3idx*60],arg7_l[41]);
    atomicAdd(&ind_arg2[42+map3idx*60],arg7_l[42]);
    atomicAdd(&ind_arg2[43+map3idx*60],arg7_l[43]);
    atomicAdd(&ind_arg2[44+map3idx*60],arg7_l[44]);
    atomicAdd(&ind_arg2[45+map3idx*60],arg7_l[45]);
    atomicAdd(&ind_arg2[46+map3idx*60],arg7_l[46]);
    atomicAdd(&ind_arg2[47+map3idx*60],arg7_l[47]);
    atomicAdd(&ind_arg2[48+map3idx*60],arg7_l[48]);
    atomicAdd(&ind_arg2[49+map3idx*60],arg7_l[49]);
    atomicAdd(&ind_arg2[50+map3idx*60],arg7_l[50]);
    atomicAdd(&ind_arg2[51+map3idx*60],arg7_l[51]);
    atomicAdd(&ind_arg2[52+map3idx*60],arg7_l[52]);
    atomicAdd(&ind_arg2[53+map3idx*60],arg7_l[53]);
    atomicAdd(&ind_arg2[54+map3idx*60],arg7_l[54]);
    atomicAdd(&ind_arg2[55+map3idx*60],arg7_l[55]);
    atomicAdd(&ind_arg2[56+map3idx*60],arg7_l[56]);
    atomicAdd(&ind_arg2[57+map3idx*60],arg7_l[57]);
    atomicAdd(&ind_arg2[58+map3idx*60],arg7_l[58]);
    atomicAdd(&ind_arg2[59+map3idx*60],arg7_l[59]);
  }
}


//host stub function
void op_par_loop_gauss_gfi_faces(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6){

  int nargs = 8;
  op_arg args[8];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 60, "double", OP_INC);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 60, "double", OP_INC);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 60, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(14);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[14].name      = name;
  OP_kernels[14].count    += 1;


  int    ninds   = 3;
  int    inds[8] = {-1,-1,0,0,1,1,2,2};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_gfi_faces\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_14
      int nthread = OP_BLOCK_SIZE_14;
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
        op_cuda_gauss_gfi_faces<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg4.data_d,
        (double *)arg6.data_d,
        arg2.map_data_d,
        (int*)arg0.data_d,
        (bool*)arg1.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[14].time     += wall_t2 - wall_t1;
}