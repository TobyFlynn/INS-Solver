//
// auto-generated by op2.py
//

//user function
__device__ void ls_flux_gpu( const int *edgeNum, const double **x, const double **y,
                    const double **sJ, const double **nx, const double **ny,
                    const double **s, double **dsldx, double **dsrdx,
                    double **dsldy, double **dsrdy) {

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
      rInd = exIndR + 7 + i - 1;
    } else {
      rInd = exIndR + i;
    }

    if(nx[0][lInd] >= 0.0) {
      dsldx[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * nx[0][lInd] * s[0][lInd];
      dsrdx[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * nx[0][lInd] * s[1][rInd];
    } else {
      dsldx[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * nx[0][lInd] * s[1][rInd];
      dsrdx[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * nx[0][lInd] * s[0][lInd];
    }

    if(ny[0][lInd] >= 0.0) {
      dsldy[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * ny[0][lInd] * s[0][lInd];
      dsrdy[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * ny[0][lInd] * s[1][rInd];
    } else {
      dsldy[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * ny[0][lInd] * s[1][rInd];
      dsrdy[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * ny[0][lInd] * s[0][lInd];
    }
  }

  for(int i = 0; i < 7; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 7 + i - 1;
    } else {
      lInd = exIndL + i;
    }

    if(nx[1][rInd] >= 0.0) {
      dsldx[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * nx[1][rInd] * s[1][rInd];
      dsrdx[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * nx[1][rInd] * s[0][lInd];
    } else {
      dsldx[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * nx[1][rInd] * s[0][lInd];
      dsrdx[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * nx[1][rInd] * s[1][rInd];
    }

    if(ny[1][rInd] >= 0.0) {
      dsldy[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * ny[1][rInd] * s[1][rInd];
      dsrdy[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * ny[1][rInd] * s[0][lInd];
    } else {
      dsldy[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * ny[1][rInd] * s[0][lInd];
      dsrdy[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * ny[1][rInd] * s[1][rInd];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_ls_flux(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  double *__restrict ind_arg6,
  double *__restrict ind_arg7,
  double *__restrict ind_arg8,
  double *__restrict ind_arg9,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  int start,
  int end,
  int   set_size) {
  double arg13_l[21];
  double arg14_l[21];
  double arg15_l[21];
  double arg16_l[21];
  double arg17_l[21];
  double arg18_l[21];
  double arg19_l[21];
  double arg20_l[21];
  double *arg13_vec[2] = {
    arg13_l,
    arg14_l,
  };
  double *arg15_vec[2] = {
    arg15_l,
    arg16_l,
  };
  double *arg17_vec[2] = {
    arg17_l,
    arg18_l,
  };
  double *arg19_vec[2] = {
    arg19_l,
    arg20_l,
  };
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg13_l[21];
    for ( int d=0; d<21; d++ ){
      arg13_l[d] = ZERO_double;
    }
    double arg14_l[21];
    for ( int d=0; d<21; d++ ){
      arg14_l[d] = ZERO_double;
    }
    double arg15_l[21];
    for ( int d=0; d<21; d++ ){
      arg15_l[d] = ZERO_double;
    }
    double arg16_l[21];
    for ( int d=0; d<21; d++ ){
      arg16_l[d] = ZERO_double;
    }
    double arg17_l[21];
    for ( int d=0; d<21; d++ ){
      arg17_l[d] = ZERO_double;
    }
    double arg18_l[21];
    for ( int d=0; d<21; d++ ){
      arg18_l[d] = ZERO_double;
    }
    double arg19_l[21];
    for ( int d=0; d<21; d++ ){
      arg19_l[d] = ZERO_double;
    }
    double arg20_l[21];
    for ( int d=0; d<21; d++ ){
      arg20_l[d] = ZERO_double;
    }
    int map1idx;
    int map2idx;
    map1idx = opDat1Map[n + set_size * 0];
    map2idx = opDat1Map[n + set_size * 1];
    const double* arg1_vec[] = {
       &ind_arg0[3 * map1idx],
       &ind_arg0[3 * map2idx]};
    const double* arg3_vec[] = {
       &ind_arg1[3 * map1idx],
       &ind_arg1[3 * map2idx]};
    const double* arg5_vec[] = {
       &ind_arg2[21 * map1idx],
       &ind_arg2[21 * map2idx]};
    const double* arg7_vec[] = {
       &ind_arg3[21 * map1idx],
       &ind_arg3[21 * map2idx]};
    const double* arg9_vec[] = {
       &ind_arg4[21 * map1idx],
       &ind_arg4[21 * map2idx]};
    const double* arg11_vec[] = {
       &ind_arg5[21 * map1idx],
       &ind_arg5[21 * map2idx]};
    double* arg13_vec[] = {
       &ind_arg6[21 * map1idx],
       &ind_arg6[21 * map2idx]};
    double* arg15_vec[] = {
       &ind_arg7[21 * map1idx],
       &ind_arg7[21 * map2idx]};
    double* arg17_vec[] = {
       &ind_arg8[21 * map1idx],
       &ind_arg8[21 * map2idx]};
    double* arg19_vec[] = {
       &ind_arg9[21 * map1idx],
       &ind_arg9[21 * map2idx]};

    //user-supplied kernel call
    ls_flux_gpu(arg0+n*2,
            arg1_vec,
            arg3_vec,
            arg5_vec,
            arg7_vec,
            arg9_vec,
            arg11_vec,
            arg13_vec,
            arg15_vec,
            arg17_vec,
            arg19_vec);
    atomicAdd(&ind_arg6[0+map1idx*21],arg13_l[0]);
    atomicAdd(&ind_arg6[1+map1idx*21],arg13_l[1]);
    atomicAdd(&ind_arg6[2+map1idx*21],arg13_l[2]);
    atomicAdd(&ind_arg6[3+map1idx*21],arg13_l[3]);
    atomicAdd(&ind_arg6[4+map1idx*21],arg13_l[4]);
    atomicAdd(&ind_arg6[5+map1idx*21],arg13_l[5]);
    atomicAdd(&ind_arg6[6+map1idx*21],arg13_l[6]);
    atomicAdd(&ind_arg6[7+map1idx*21],arg13_l[7]);
    atomicAdd(&ind_arg6[8+map1idx*21],arg13_l[8]);
    atomicAdd(&ind_arg6[9+map1idx*21],arg13_l[9]);
    atomicAdd(&ind_arg6[10+map1idx*21],arg13_l[10]);
    atomicAdd(&ind_arg6[11+map1idx*21],arg13_l[11]);
    atomicAdd(&ind_arg6[12+map1idx*21],arg13_l[12]);
    atomicAdd(&ind_arg6[13+map1idx*21],arg13_l[13]);
    atomicAdd(&ind_arg6[14+map1idx*21],arg13_l[14]);
    atomicAdd(&ind_arg6[15+map1idx*21],arg13_l[15]);
    atomicAdd(&ind_arg6[16+map1idx*21],arg13_l[16]);
    atomicAdd(&ind_arg6[17+map1idx*21],arg13_l[17]);
    atomicAdd(&ind_arg6[18+map1idx*21],arg13_l[18]);
    atomicAdd(&ind_arg6[19+map1idx*21],arg13_l[19]);
    atomicAdd(&ind_arg6[20+map1idx*21],arg13_l[20]);
    atomicAdd(&ind_arg6[0+map2idx*21],arg14_l[0]);
    atomicAdd(&ind_arg6[1+map2idx*21],arg14_l[1]);
    atomicAdd(&ind_arg6[2+map2idx*21],arg14_l[2]);
    atomicAdd(&ind_arg6[3+map2idx*21],arg14_l[3]);
    atomicAdd(&ind_arg6[4+map2idx*21],arg14_l[4]);
    atomicAdd(&ind_arg6[5+map2idx*21],arg14_l[5]);
    atomicAdd(&ind_arg6[6+map2idx*21],arg14_l[6]);
    atomicAdd(&ind_arg6[7+map2idx*21],arg14_l[7]);
    atomicAdd(&ind_arg6[8+map2idx*21],arg14_l[8]);
    atomicAdd(&ind_arg6[9+map2idx*21],arg14_l[9]);
    atomicAdd(&ind_arg6[10+map2idx*21],arg14_l[10]);
    atomicAdd(&ind_arg6[11+map2idx*21],arg14_l[11]);
    atomicAdd(&ind_arg6[12+map2idx*21],arg14_l[12]);
    atomicAdd(&ind_arg6[13+map2idx*21],arg14_l[13]);
    atomicAdd(&ind_arg6[14+map2idx*21],arg14_l[14]);
    atomicAdd(&ind_arg6[15+map2idx*21],arg14_l[15]);
    atomicAdd(&ind_arg6[16+map2idx*21],arg14_l[16]);
    atomicAdd(&ind_arg6[17+map2idx*21],arg14_l[17]);
    atomicAdd(&ind_arg6[18+map2idx*21],arg14_l[18]);
    atomicAdd(&ind_arg6[19+map2idx*21],arg14_l[19]);
    atomicAdd(&ind_arg6[20+map2idx*21],arg14_l[20]);
    atomicAdd(&ind_arg7[0+map1idx*21],arg15_l[0]);
    atomicAdd(&ind_arg7[1+map1idx*21],arg15_l[1]);
    atomicAdd(&ind_arg7[2+map1idx*21],arg15_l[2]);
    atomicAdd(&ind_arg7[3+map1idx*21],arg15_l[3]);
    atomicAdd(&ind_arg7[4+map1idx*21],arg15_l[4]);
    atomicAdd(&ind_arg7[5+map1idx*21],arg15_l[5]);
    atomicAdd(&ind_arg7[6+map1idx*21],arg15_l[6]);
    atomicAdd(&ind_arg7[7+map1idx*21],arg15_l[7]);
    atomicAdd(&ind_arg7[8+map1idx*21],arg15_l[8]);
    atomicAdd(&ind_arg7[9+map1idx*21],arg15_l[9]);
    atomicAdd(&ind_arg7[10+map1idx*21],arg15_l[10]);
    atomicAdd(&ind_arg7[11+map1idx*21],arg15_l[11]);
    atomicAdd(&ind_arg7[12+map1idx*21],arg15_l[12]);
    atomicAdd(&ind_arg7[13+map1idx*21],arg15_l[13]);
    atomicAdd(&ind_arg7[14+map1idx*21],arg15_l[14]);
    atomicAdd(&ind_arg7[15+map1idx*21],arg15_l[15]);
    atomicAdd(&ind_arg7[16+map1idx*21],arg15_l[16]);
    atomicAdd(&ind_arg7[17+map1idx*21],arg15_l[17]);
    atomicAdd(&ind_arg7[18+map1idx*21],arg15_l[18]);
    atomicAdd(&ind_arg7[19+map1idx*21],arg15_l[19]);
    atomicAdd(&ind_arg7[20+map1idx*21],arg15_l[20]);
    atomicAdd(&ind_arg7[0+map2idx*21],arg16_l[0]);
    atomicAdd(&ind_arg7[1+map2idx*21],arg16_l[1]);
    atomicAdd(&ind_arg7[2+map2idx*21],arg16_l[2]);
    atomicAdd(&ind_arg7[3+map2idx*21],arg16_l[3]);
    atomicAdd(&ind_arg7[4+map2idx*21],arg16_l[4]);
    atomicAdd(&ind_arg7[5+map2idx*21],arg16_l[5]);
    atomicAdd(&ind_arg7[6+map2idx*21],arg16_l[6]);
    atomicAdd(&ind_arg7[7+map2idx*21],arg16_l[7]);
    atomicAdd(&ind_arg7[8+map2idx*21],arg16_l[8]);
    atomicAdd(&ind_arg7[9+map2idx*21],arg16_l[9]);
    atomicAdd(&ind_arg7[10+map2idx*21],arg16_l[10]);
    atomicAdd(&ind_arg7[11+map2idx*21],arg16_l[11]);
    atomicAdd(&ind_arg7[12+map2idx*21],arg16_l[12]);
    atomicAdd(&ind_arg7[13+map2idx*21],arg16_l[13]);
    atomicAdd(&ind_arg7[14+map2idx*21],arg16_l[14]);
    atomicAdd(&ind_arg7[15+map2idx*21],arg16_l[15]);
    atomicAdd(&ind_arg7[16+map2idx*21],arg16_l[16]);
    atomicAdd(&ind_arg7[17+map2idx*21],arg16_l[17]);
    atomicAdd(&ind_arg7[18+map2idx*21],arg16_l[18]);
    atomicAdd(&ind_arg7[19+map2idx*21],arg16_l[19]);
    atomicAdd(&ind_arg7[20+map2idx*21],arg16_l[20]);
    atomicAdd(&ind_arg8[0+map1idx*21],arg17_l[0]);
    atomicAdd(&ind_arg8[1+map1idx*21],arg17_l[1]);
    atomicAdd(&ind_arg8[2+map1idx*21],arg17_l[2]);
    atomicAdd(&ind_arg8[3+map1idx*21],arg17_l[3]);
    atomicAdd(&ind_arg8[4+map1idx*21],arg17_l[4]);
    atomicAdd(&ind_arg8[5+map1idx*21],arg17_l[5]);
    atomicAdd(&ind_arg8[6+map1idx*21],arg17_l[6]);
    atomicAdd(&ind_arg8[7+map1idx*21],arg17_l[7]);
    atomicAdd(&ind_arg8[8+map1idx*21],arg17_l[8]);
    atomicAdd(&ind_arg8[9+map1idx*21],arg17_l[9]);
    atomicAdd(&ind_arg8[10+map1idx*21],arg17_l[10]);
    atomicAdd(&ind_arg8[11+map1idx*21],arg17_l[11]);
    atomicAdd(&ind_arg8[12+map1idx*21],arg17_l[12]);
    atomicAdd(&ind_arg8[13+map1idx*21],arg17_l[13]);
    atomicAdd(&ind_arg8[14+map1idx*21],arg17_l[14]);
    atomicAdd(&ind_arg8[15+map1idx*21],arg17_l[15]);
    atomicAdd(&ind_arg8[16+map1idx*21],arg17_l[16]);
    atomicAdd(&ind_arg8[17+map1idx*21],arg17_l[17]);
    atomicAdd(&ind_arg8[18+map1idx*21],arg17_l[18]);
    atomicAdd(&ind_arg8[19+map1idx*21],arg17_l[19]);
    atomicAdd(&ind_arg8[20+map1idx*21],arg17_l[20]);
    atomicAdd(&ind_arg8[0+map2idx*21],arg18_l[0]);
    atomicAdd(&ind_arg8[1+map2idx*21],arg18_l[1]);
    atomicAdd(&ind_arg8[2+map2idx*21],arg18_l[2]);
    atomicAdd(&ind_arg8[3+map2idx*21],arg18_l[3]);
    atomicAdd(&ind_arg8[4+map2idx*21],arg18_l[4]);
    atomicAdd(&ind_arg8[5+map2idx*21],arg18_l[5]);
    atomicAdd(&ind_arg8[6+map2idx*21],arg18_l[6]);
    atomicAdd(&ind_arg8[7+map2idx*21],arg18_l[7]);
    atomicAdd(&ind_arg8[8+map2idx*21],arg18_l[8]);
    atomicAdd(&ind_arg8[9+map2idx*21],arg18_l[9]);
    atomicAdd(&ind_arg8[10+map2idx*21],arg18_l[10]);
    atomicAdd(&ind_arg8[11+map2idx*21],arg18_l[11]);
    atomicAdd(&ind_arg8[12+map2idx*21],arg18_l[12]);
    atomicAdd(&ind_arg8[13+map2idx*21],arg18_l[13]);
    atomicAdd(&ind_arg8[14+map2idx*21],arg18_l[14]);
    atomicAdd(&ind_arg8[15+map2idx*21],arg18_l[15]);
    atomicAdd(&ind_arg8[16+map2idx*21],arg18_l[16]);
    atomicAdd(&ind_arg8[17+map2idx*21],arg18_l[17]);
    atomicAdd(&ind_arg8[18+map2idx*21],arg18_l[18]);
    atomicAdd(&ind_arg8[19+map2idx*21],arg18_l[19]);
    atomicAdd(&ind_arg8[20+map2idx*21],arg18_l[20]);
    atomicAdd(&ind_arg9[0+map1idx*21],arg19_l[0]);
    atomicAdd(&ind_arg9[1+map1idx*21],arg19_l[1]);
    atomicAdd(&ind_arg9[2+map1idx*21],arg19_l[2]);
    atomicAdd(&ind_arg9[3+map1idx*21],arg19_l[3]);
    atomicAdd(&ind_arg9[4+map1idx*21],arg19_l[4]);
    atomicAdd(&ind_arg9[5+map1idx*21],arg19_l[5]);
    atomicAdd(&ind_arg9[6+map1idx*21],arg19_l[6]);
    atomicAdd(&ind_arg9[7+map1idx*21],arg19_l[7]);
    atomicAdd(&ind_arg9[8+map1idx*21],arg19_l[8]);
    atomicAdd(&ind_arg9[9+map1idx*21],arg19_l[9]);
    atomicAdd(&ind_arg9[10+map1idx*21],arg19_l[10]);
    atomicAdd(&ind_arg9[11+map1idx*21],arg19_l[11]);
    atomicAdd(&ind_arg9[12+map1idx*21],arg19_l[12]);
    atomicAdd(&ind_arg9[13+map1idx*21],arg19_l[13]);
    atomicAdd(&ind_arg9[14+map1idx*21],arg19_l[14]);
    atomicAdd(&ind_arg9[15+map1idx*21],arg19_l[15]);
    atomicAdd(&ind_arg9[16+map1idx*21],arg19_l[16]);
    atomicAdd(&ind_arg9[17+map1idx*21],arg19_l[17]);
    atomicAdd(&ind_arg9[18+map1idx*21],arg19_l[18]);
    atomicAdd(&ind_arg9[19+map1idx*21],arg19_l[19]);
    atomicAdd(&ind_arg9[20+map1idx*21],arg19_l[20]);
    atomicAdd(&ind_arg9[0+map2idx*21],arg20_l[0]);
    atomicAdd(&ind_arg9[1+map2idx*21],arg20_l[1]);
    atomicAdd(&ind_arg9[2+map2idx*21],arg20_l[2]);
    atomicAdd(&ind_arg9[3+map2idx*21],arg20_l[3]);
    atomicAdd(&ind_arg9[4+map2idx*21],arg20_l[4]);
    atomicAdd(&ind_arg9[5+map2idx*21],arg20_l[5]);
    atomicAdd(&ind_arg9[6+map2idx*21],arg20_l[6]);
    atomicAdd(&ind_arg9[7+map2idx*21],arg20_l[7]);
    atomicAdd(&ind_arg9[8+map2idx*21],arg20_l[8]);
    atomicAdd(&ind_arg9[9+map2idx*21],arg20_l[9]);
    atomicAdd(&ind_arg9[10+map2idx*21],arg20_l[10]);
    atomicAdd(&ind_arg9[11+map2idx*21],arg20_l[11]);
    atomicAdd(&ind_arg9[12+map2idx*21],arg20_l[12]);
    atomicAdd(&ind_arg9[13+map2idx*21],arg20_l[13]);
    atomicAdd(&ind_arg9[14+map2idx*21],arg20_l[14]);
    atomicAdd(&ind_arg9[15+map2idx*21],arg20_l[15]);
    atomicAdd(&ind_arg9[16+map2idx*21],arg20_l[16]);
    atomicAdd(&ind_arg9[17+map2idx*21],arg20_l[17]);
    atomicAdd(&ind_arg9[18+map2idx*21],arg20_l[18]);
    atomicAdd(&ind_arg9[19+map2idx*21],arg20_l[19]);
    atomicAdd(&ind_arg9[20+map2idx*21],arg20_l[20]);
  }
}


//host stub function
void op_par_loop_ls_flux(char const *name, op_set set,
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
  op_arg arg19){

  int nargs = 21;
  op_arg args[21];

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

  arg17.idx = 0;
  args[17] = arg17;
  for ( int v=1; v<2; v++ ){
    args[17 + v] = op_arg_dat(arg17.dat, v, arg17.map, 21, "double", OP_INC);
  }

  arg19.idx = 0;
  args[19] = arg19;
  for ( int v=1; v<2; v++ ){
    args[19 + v] = op_arg_dat(arg19.dat, v, arg19.map, 21, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(57);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[57].name      = name;
  OP_kernels[57].count    += 1;


  int    ninds   = 10;
  int    inds[21] = {-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_flux\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_57
      int nthread = OP_BLOCK_SIZE_57;
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
        op_cuda_ls_flux<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg3.data_d,
        (double *)arg5.data_d,
        (double *)arg7.data_d,
        (double *)arg9.data_d,
        (double *)arg11.data_d,
        (double *)arg13.data_d,
        (double *)arg15.data_d,
        (double *)arg17.data_d,
        (double *)arg19.data_d,
        arg1.map_data_d,
        (int*)arg0.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[57].time     += wall_t2 - wall_t1;
}
