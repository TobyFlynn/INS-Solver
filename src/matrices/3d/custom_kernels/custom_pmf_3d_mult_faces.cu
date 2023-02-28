/*
template<int dg_npf>
__device__ void _pmf_3d_mult_faces_gpu(const int *faceNum, const int *fmaskL_corrected,
                              const int *fmaskR_corrected, const DG_FP *nx,
                              const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                              const DG_FP *sJ, const DG_FP **in, const DG_FP **in_x,
                              const DG_FP **in_y, const DG_FP **in_z, DG_FP **l_x,
                              DG_FP **l_y, DG_FP **l_z, DG_FP **out) {
  const int findL = faceNum[0] * dg_npf;
  const int findR = faceNum[1] * dg_npf;
  const int FMASK__[][] = {{0, 1, 2, 0, 1, 3, 1, 2, 3, 0, 2, 3},
                            {0, 1, 2, 3, 4, 5, 0, 1, 2, 6, 7, 9, 2, 4, 5, 7, 8, 9, 0, 3, 5, 6, 8, 9},
                            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 10, 11, 12, 16, 17, 19, 3, 6, 8, 9, 12, 14, 15, 17, 18, 19, 0, 4, 7, 9, 10, 13, 15, 16, 18, 19}};
  const int *fmask = FMASK__[p];
  const int *fmaskL = &fmask[faceNum[0] * dg_npf];
  const int *fmaskR = &fmask[faceNum[1] * dg_npf];

  const DG_FP gtau = 2.0 * (3 + 1) * (3 + 1) * fmax(fscale[0], fscale[1]);

  for(int j = 0; j < dg_npf; j++) {
    const DG_FP diffL_u = in[0][fmaskL[j]] - in[1][fmaskR_corrected[j]];
    const DG_FP diffL_u_x = nx[0] * (in_x[1][fmaskR_corrected[j]] + in_x[0][fmaskL[j]]);
    const DG_FP diffL_u_y = ny[0] * (in_y[1][fmaskR_corrected[j]] + in_y[0][fmaskL[j]]);
    const DG_FP diffL_u_z = nz[0] * (in_z[1][fmaskR_corrected[j]] + in_z[0][fmaskL[j]]);
    const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

    const int indL = findL + j;
    out[0][indL] += 0.5 * sJ[0] * (gtau * diffL_u - diffL_u_grad);
    const DG_FP l_tmpL = 0.5 * sJ[0] * -diffL_u;
    l_x[0][indL] += nx[0] * l_tmpL;
    l_y[0][indL] += ny[0] * l_tmpL;
    l_z[0][indL] += nz[0] * l_tmpL;

    const DG_FP diffR_u = in[1][fmaskR[j]] - in[0][fmaskL_corrected[j]];
    const DG_FP diffR_u_x = nx[1] * (in_x[1][fmaskR[j]] + in_x[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_y = ny[1] * (in_y[1][fmaskR[j]] + in_y[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_z = nz[1] * (in_z[1][fmaskR[j]] + in_z[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_grad = diffR_u_x + diffR_u_y + diffR_u_z;

    const int indR = findR + j;
    out[1][indR] += 0.5 * sJ[1] * (gtau * diffR_u - diffR_u_grad);
    const DG_FP l_tmpR = 0.5 * sJ[1] * -diffR_u;
    l_x[1][indR] += nx[1] * l_tmpR;
    l_y[1][indR] += ny[1] * l_tmpR;
    l_z[1][indR] += nz[1] * l_tmpR;
  }

}
*/
__device__ void _1_pmf_3d_mult_faces_gpu(const int *faceNum, const int *fmaskL_corrected,
                              const int *fmaskR_corrected, const DG_FP *nx,
                              const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                              const DG_FP *sJ, const DG_FP **in, const DG_FP **in_x,
                              const DG_FP **in_y, const DG_FP **in_z, DG_FP **l_x,
                              DG_FP **l_y, DG_FP **l_z, DG_FP **out) {
  const int findL = faceNum[0] * 3;
  const int findR = faceNum[1] * 3;
  const int fmask[] = {0, 1, 2, 0, 1, 3, 1, 2, 3, 0, 2, 3};
  const int *fmaskL = &fmask[faceNum[0] * 3];
  const int *fmaskR = &fmask[faceNum[1] * 3];

  const DG_FP gtau = 2.0 * (3 + 1) * (3 + 1) * fmax(fscale[0], fscale[1]);

  for(int j = 0; j < 3; j++) {
    const DG_FP diffL_u = in[0][fmaskL[j]] - in[1][fmaskR_corrected[j]];
    const DG_FP diffL_u_x = nx[0] * (in_x[1][fmaskR_corrected[j]] + in_x[0][fmaskL[j]]);
    const DG_FP diffL_u_y = ny[0] * (in_y[1][fmaskR_corrected[j]] + in_y[0][fmaskL[j]]);
    const DG_FP diffL_u_z = nz[0] * (in_z[1][fmaskR_corrected[j]] + in_z[0][fmaskL[j]]);
    const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

    const int indL = findL + j;
    out[0][indL] += 0.5 * sJ[0] * (gtau * diffL_u - diffL_u_grad);
    const DG_FP l_tmpL = 0.5 * sJ[0] * -diffL_u;
    l_x[0][indL] += nx[0] * l_tmpL;
    l_y[0][indL] += ny[0] * l_tmpL;
    l_z[0][indL] += nz[0] * l_tmpL;

    const DG_FP diffR_u = in[1][fmaskR[j]] - in[0][fmaskL_corrected[j]];
    const DG_FP diffR_u_x = nx[1] * (in_x[1][fmaskR[j]] + in_x[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_y = ny[1] * (in_y[1][fmaskR[j]] + in_y[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_z = nz[1] * (in_z[1][fmaskR[j]] + in_z[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_grad = diffR_u_x + diffR_u_y + diffR_u_z;

    const int indR = findR + j;
    out[1][indR] += 0.5 * sJ[1] * (gtau * diffR_u - diffR_u_grad);
    const DG_FP l_tmpR = 0.5 * sJ[1] * -diffR_u;
    l_x[1][indR] += nx[1] * l_tmpR;
    l_y[1][indR] += ny[1] * l_tmpR;
    l_z[1][indR] += nz[1] * l_tmpR;
  }

}

__device__ void _2_pmf_3d_mult_faces_gpu(const int *faceNum, const int *fmaskL_corrected,
                              const int *fmaskR_corrected, const DG_FP *nx,
                              const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                              const DG_FP *sJ, const DG_FP **in, const DG_FP **in_x,
                              const DG_FP **in_y, const DG_FP **in_z, DG_FP **l_x,
                              DG_FP **l_y, DG_FP **l_z, DG_FP **out) {
  const int findL = faceNum[0] * 6;
  const int findR = faceNum[1] * 6;
  const int fmask[] = {0, 1, 2, 3, 4, 5, 0, 1, 2, 6, 7, 9, 2, 4, 5, 7, 8, 9, 0, 3, 5, 6, 8, 9};
  const int *fmaskL = &fmask[faceNum[0] * 6];
  const int *fmaskR = &fmask[faceNum[1] * 6];

  const DG_FP gtau = 2.0 * (3 + 1) * (3 + 1) * fmax(fscale[0], fscale[1]);

  for(int j = 0; j < 6; j++) {
    const DG_FP diffL_u = in[0][fmaskL[j]] - in[1][fmaskR_corrected[j]];
    const DG_FP diffL_u_x = nx[0] * (in_x[1][fmaskR_corrected[j]] + in_x[0][fmaskL[j]]);
    const DG_FP diffL_u_y = ny[0] * (in_y[1][fmaskR_corrected[j]] + in_y[0][fmaskL[j]]);
    const DG_FP diffL_u_z = nz[0] * (in_z[1][fmaskR_corrected[j]] + in_z[0][fmaskL[j]]);
    const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

    const int indL = findL + j;
    out[0][indL] += 0.5 * sJ[0] * (gtau * diffL_u - diffL_u_grad);
    const DG_FP l_tmpL = 0.5 * sJ[0] * -diffL_u;
    l_x[0][indL] += nx[0] * l_tmpL;
    l_y[0][indL] += ny[0] * l_tmpL;
    l_z[0][indL] += nz[0] * l_tmpL;

    const DG_FP diffR_u = in[1][fmaskR[j]] - in[0][fmaskL_corrected[j]];
    const DG_FP diffR_u_x = nx[1] * (in_x[1][fmaskR[j]] + in_x[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_y = ny[1] * (in_y[1][fmaskR[j]] + in_y[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_z = nz[1] * (in_z[1][fmaskR[j]] + in_z[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_grad = diffR_u_x + diffR_u_y + diffR_u_z;

    const int indR = findR + j;
    out[1][indR] += 0.5 * sJ[1] * (gtau * diffR_u - diffR_u_grad);
    const DG_FP l_tmpR = 0.5 * sJ[1] * -diffR_u;
    l_x[1][indR] += nx[1] * l_tmpR;
    l_y[1][indR] += ny[1] * l_tmpR;
    l_z[1][indR] += nz[1] * l_tmpR;
  }

}

__device__ void _3_pmf_3d_mult_faces_gpu(const int *faceNum, const int *fmaskL_corrected,
                              const int *fmaskR_corrected, const DG_FP *nx,
                              const DG_FP *ny, const DG_FP *nz, const DG_FP *fscale,
                              const DG_FP *sJ, const DG_FP **in, const DG_FP **in_x,
                              const DG_FP **in_y, const DG_FP **in_z, DG_FP **l_x,
                              DG_FP **l_y, DG_FP **l_z, DG_FP **out) {
  const int findL = faceNum[0] * 10;
  const int findR = faceNum[1] * 10;
  const int fmask[]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 10, 11, 12, 16, 17, 19, 3, 6, 8, 9, 12, 14, 15, 17, 18, 19, 0, 4, 7, 9, 10, 13, 15, 16, 18, 19};
  const int *fmaskL = &fmask[faceNum[0] * 10];
  const int *fmaskR = &fmask[faceNum[1] * 10];

  const DG_FP gtau = 2.0 * (3 + 1) * (3 + 1) * fmax(fscale[0], fscale[1]);

  for(int j = 0; j < 10; j++) {
    const DG_FP diffL_u = in[0][fmaskL[j]] - in[1][fmaskR_corrected[j]];
    const DG_FP diffL_u_x = nx[0] * (in_x[1][fmaskR_corrected[j]] + in_x[0][fmaskL[j]]);
    const DG_FP diffL_u_y = ny[0] * (in_y[1][fmaskR_corrected[j]] + in_y[0][fmaskL[j]]);
    const DG_FP diffL_u_z = nz[0] * (in_z[1][fmaskR_corrected[j]] + in_z[0][fmaskL[j]]);
    const DG_FP diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

    const int indL = findL + j;
    out[0][indL] += 0.5 * sJ[0] * (gtau * diffL_u - diffL_u_grad);
    const DG_FP l_tmpL = 0.5 * sJ[0] * -diffL_u;
    l_x[0][indL] += nx[0] * l_tmpL;
    l_y[0][indL] += ny[0] * l_tmpL;
    l_z[0][indL] += nz[0] * l_tmpL;

    const DG_FP diffR_u = in[1][fmaskR[j]] - in[0][fmaskL_corrected[j]];
    const DG_FP diffR_u_x = nx[1] * (in_x[1][fmaskR[j]] + in_x[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_y = ny[1] * (in_y[1][fmaskR[j]] + in_y[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_z = nz[1] * (in_z[1][fmaskR[j]] + in_z[0][fmaskL_corrected[j]]);
    const DG_FP diffR_u_grad = diffR_u_x + diffR_u_y + diffR_u_z;

    const int indR = findR + j;
    out[1][indR] += 0.5 * sJ[1] * (gtau * diffR_u - diffR_u_grad);
    const DG_FP l_tmpR = 0.5 * sJ[1] * -diffR_u;
    l_x[1][indR] += nx[1] * l_tmpR;
    l_y[1][indR] += ny[1] * l_tmpR;
    l_z[1][indR] += nz[1] * l_tmpR;
  }

}

// CUDA kernel function
__global__ void _op_cuda_pmf_3d_mult_faces(
  const int *__restrict ind_arg0,
  const DG_FP *__restrict ind_arg1,
  const DG_FP *__restrict ind_arg2,
  const DG_FP *__restrict ind_arg3,
  const DG_FP *__restrict ind_arg4,
  DG_FP *__restrict ind_arg5,
  DG_FP *__restrict ind_arg6,
  DG_FP *__restrict ind_arg7,
  DG_FP *__restrict ind_arg8,
  const int *__restrict opDat0Map,
  const int *__restrict arg2,
  const int *__restrict arg3,
  const int *__restrict arg4,
  const DG_FP *__restrict arg5,
  const DG_FP *__restrict arg6,
  const DG_FP *__restrict arg7,
  const DG_FP *__restrict arg8,
  const DG_FP *__restrict arg9,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg18_l[40];
    for ( int d=0; d<40; d++ ){
      arg18_l[d] = ZERO_double;
    }
    double arg19_l[40];
    for ( int d=0; d<40; d++ ){
      arg19_l[d] = ZERO_double;
    }
    double arg20_l[40];
    for ( int d=0; d<40; d++ ){
      arg20_l[d] = ZERO_double;
    }
    double arg21_l[40];
    for ( int d=0; d<40; d++ ){
      arg21_l[d] = ZERO_double;
    }
    double arg22_l[40];
    for ( int d=0; d<40; d++ ){
      arg22_l[d] = ZERO_double;
    }
    double arg23_l[40];
    for ( int d=0; d<40; d++ ){
      arg23_l[d] = ZERO_double;
    }
    double arg24_l[40];
    for ( int d=0; d<40; d++ ){
      arg24_l[d] = ZERO_double;
    }
    double arg25_l[40];
    for ( int d=0; d<40; d++ ){
      arg25_l[d] = ZERO_double;
    }
    int map0idx;
    int map1idx;
    map0idx = opDat0Map[n + set_size * 0];
    map1idx = opDat0Map[n + set_size * 1];
    const int* arg0_vec[] = {
       &ind_arg0[1 * map0idx],
       &ind_arg0[1 * map1idx]};
    const double* arg10_vec[] = {
       &ind_arg1[20 * map0idx],
       &ind_arg1[20 * map1idx]};
    const double* arg12_vec[] = {
       &ind_arg2[20 * map0idx],
       &ind_arg2[20 * map1idx]};
    const double* arg14_vec[] = {
       &ind_arg3[20 * map0idx],
       &ind_arg3[20 * map1idx]};
    const double* arg16_vec[] = {
       &ind_arg4[20 * map0idx],
       &ind_arg4[20 * map1idx]};
    double* arg18_vec[] = {
      arg18_l,
      arg19_l};
    double* arg20_vec[] = {
      arg20_l,
      arg21_l};
    double* arg22_vec[] = {
      arg22_l,
      arg23_l};
    double* arg24_vec[] = {
      arg24_l,
      arg25_l};

    //user-supplied kernel call
    switch(arg0_vec[0][0]) {
      case 1:
        _1_pmf_3d_mult_faces_gpu(
                          arg2+n*2,
                          arg3+n*10,
                          arg4+n*10,
                          arg5+n*2,
                          arg6+n*2,
                          arg7+n*2,
                          arg8+n*2,
                          arg9+n*2,
                          arg10_vec,
                          arg12_vec,
                          arg14_vec,
                          arg16_vec,
                          arg18_vec,
                          arg20_vec,
                          arg22_vec,
                          arg24_vec);
        break;
      case 2:
        _2_pmf_3d_mult_faces_gpu(
                          arg2+n*2,
                          arg3+n*10,
                          arg4+n*10,
                          arg5+n*2,
                          arg6+n*2,
                          arg7+n*2,
                          arg8+n*2,
                          arg9+n*2,
                          arg10_vec,
                          arg12_vec,
                          arg14_vec,
                          arg16_vec,
                          arg18_vec,
                          arg20_vec,
                          arg22_vec,
                          arg24_vec);
        break;
      case 3:
        _3_pmf_3d_mult_faces_gpu(
                          arg2+n*2,
                          arg3+n*10,
                          arg4+n*10,
                          arg5+n*2,
                          arg6+n*2,
                          arg7+n*2,
                          arg8+n*2,
                          arg9+n*2,
                          arg10_vec,
                          arg12_vec,
                          arg14_vec,
                          arg16_vec,
                          arg18_vec,
                          arg20_vec,
                          arg22_vec,
                          arg24_vec);
        break;
    }

    atomicAdd(&ind_arg5[0+map0idx*40],arg18_l[0]);
    atomicAdd(&ind_arg5[1+map0idx*40],arg18_l[1]);
    atomicAdd(&ind_arg5[2+map0idx*40],arg18_l[2]);
    atomicAdd(&ind_arg5[3+map0idx*40],arg18_l[3]);
    atomicAdd(&ind_arg5[4+map0idx*40],arg18_l[4]);
    atomicAdd(&ind_arg5[5+map0idx*40],arg18_l[5]);
    atomicAdd(&ind_arg5[6+map0idx*40],arg18_l[6]);
    atomicAdd(&ind_arg5[7+map0idx*40],arg18_l[7]);
    atomicAdd(&ind_arg5[8+map0idx*40],arg18_l[8]);
    atomicAdd(&ind_arg5[9+map0idx*40],arg18_l[9]);
    atomicAdd(&ind_arg5[10+map0idx*40],arg18_l[10]);
    atomicAdd(&ind_arg5[11+map0idx*40],arg18_l[11]);
    atomicAdd(&ind_arg5[12+map0idx*40],arg18_l[12]);
    atomicAdd(&ind_arg5[13+map0idx*40],arg18_l[13]);
    atomicAdd(&ind_arg5[14+map0idx*40],arg18_l[14]);
    atomicAdd(&ind_arg5[15+map0idx*40],arg18_l[15]);
    atomicAdd(&ind_arg5[16+map0idx*40],arg18_l[16]);
    atomicAdd(&ind_arg5[17+map0idx*40],arg18_l[17]);
    atomicAdd(&ind_arg5[18+map0idx*40],arg18_l[18]);
    atomicAdd(&ind_arg5[19+map0idx*40],arg18_l[19]);
    atomicAdd(&ind_arg5[20+map0idx*40],arg18_l[20]);
    atomicAdd(&ind_arg5[21+map0idx*40],arg18_l[21]);
    atomicAdd(&ind_arg5[22+map0idx*40],arg18_l[22]);
    atomicAdd(&ind_arg5[23+map0idx*40],arg18_l[23]);
    atomicAdd(&ind_arg5[24+map0idx*40],arg18_l[24]);
    atomicAdd(&ind_arg5[25+map0idx*40],arg18_l[25]);
    atomicAdd(&ind_arg5[26+map0idx*40],arg18_l[26]);
    atomicAdd(&ind_arg5[27+map0idx*40],arg18_l[27]);
    atomicAdd(&ind_arg5[28+map0idx*40],arg18_l[28]);
    atomicAdd(&ind_arg5[29+map0idx*40],arg18_l[29]);
    atomicAdd(&ind_arg5[30+map0idx*40],arg18_l[30]);
    atomicAdd(&ind_arg5[31+map0idx*40],arg18_l[31]);
    atomicAdd(&ind_arg5[32+map0idx*40],arg18_l[32]);
    atomicAdd(&ind_arg5[33+map0idx*40],arg18_l[33]);
    atomicAdd(&ind_arg5[34+map0idx*40],arg18_l[34]);
    atomicAdd(&ind_arg5[35+map0idx*40],arg18_l[35]);
    atomicAdd(&ind_arg5[36+map0idx*40],arg18_l[36]);
    atomicAdd(&ind_arg5[37+map0idx*40],arg18_l[37]);
    atomicAdd(&ind_arg5[38+map0idx*40],arg18_l[38]);
    atomicAdd(&ind_arg5[39+map0idx*40],arg18_l[39]);
    atomicAdd(&ind_arg5[0+map1idx*40],arg19_l[0]);
    atomicAdd(&ind_arg5[1+map1idx*40],arg19_l[1]);
    atomicAdd(&ind_arg5[2+map1idx*40],arg19_l[2]);
    atomicAdd(&ind_arg5[3+map1idx*40],arg19_l[3]);
    atomicAdd(&ind_arg5[4+map1idx*40],arg19_l[4]);
    atomicAdd(&ind_arg5[5+map1idx*40],arg19_l[5]);
    atomicAdd(&ind_arg5[6+map1idx*40],arg19_l[6]);
    atomicAdd(&ind_arg5[7+map1idx*40],arg19_l[7]);
    atomicAdd(&ind_arg5[8+map1idx*40],arg19_l[8]);
    atomicAdd(&ind_arg5[9+map1idx*40],arg19_l[9]);
    atomicAdd(&ind_arg5[10+map1idx*40],arg19_l[10]);
    atomicAdd(&ind_arg5[11+map1idx*40],arg19_l[11]);
    atomicAdd(&ind_arg5[12+map1idx*40],arg19_l[12]);
    atomicAdd(&ind_arg5[13+map1idx*40],arg19_l[13]);
    atomicAdd(&ind_arg5[14+map1idx*40],arg19_l[14]);
    atomicAdd(&ind_arg5[15+map1idx*40],arg19_l[15]);
    atomicAdd(&ind_arg5[16+map1idx*40],arg19_l[16]);
    atomicAdd(&ind_arg5[17+map1idx*40],arg19_l[17]);
    atomicAdd(&ind_arg5[18+map1idx*40],arg19_l[18]);
    atomicAdd(&ind_arg5[19+map1idx*40],arg19_l[19]);
    atomicAdd(&ind_arg5[20+map1idx*40],arg19_l[20]);
    atomicAdd(&ind_arg5[21+map1idx*40],arg19_l[21]);
    atomicAdd(&ind_arg5[22+map1idx*40],arg19_l[22]);
    atomicAdd(&ind_arg5[23+map1idx*40],arg19_l[23]);
    atomicAdd(&ind_arg5[24+map1idx*40],arg19_l[24]);
    atomicAdd(&ind_arg5[25+map1idx*40],arg19_l[25]);
    atomicAdd(&ind_arg5[26+map1idx*40],arg19_l[26]);
    atomicAdd(&ind_arg5[27+map1idx*40],arg19_l[27]);
    atomicAdd(&ind_arg5[28+map1idx*40],arg19_l[28]);
    atomicAdd(&ind_arg5[29+map1idx*40],arg19_l[29]);
    atomicAdd(&ind_arg5[30+map1idx*40],arg19_l[30]);
    atomicAdd(&ind_arg5[31+map1idx*40],arg19_l[31]);
    atomicAdd(&ind_arg5[32+map1idx*40],arg19_l[32]);
    atomicAdd(&ind_arg5[33+map1idx*40],arg19_l[33]);
    atomicAdd(&ind_arg5[34+map1idx*40],arg19_l[34]);
    atomicAdd(&ind_arg5[35+map1idx*40],arg19_l[35]);
    atomicAdd(&ind_arg5[36+map1idx*40],arg19_l[36]);
    atomicAdd(&ind_arg5[37+map1idx*40],arg19_l[37]);
    atomicAdd(&ind_arg5[38+map1idx*40],arg19_l[38]);
    atomicAdd(&ind_arg5[39+map1idx*40],arg19_l[39]);
    atomicAdd(&ind_arg6[0+map0idx*40],arg20_l[0]);
    atomicAdd(&ind_arg6[1+map0idx*40],arg20_l[1]);
    atomicAdd(&ind_arg6[2+map0idx*40],arg20_l[2]);
    atomicAdd(&ind_arg6[3+map0idx*40],arg20_l[3]);
    atomicAdd(&ind_arg6[4+map0idx*40],arg20_l[4]);
    atomicAdd(&ind_arg6[5+map0idx*40],arg20_l[5]);
    atomicAdd(&ind_arg6[6+map0idx*40],arg20_l[6]);
    atomicAdd(&ind_arg6[7+map0idx*40],arg20_l[7]);
    atomicAdd(&ind_arg6[8+map0idx*40],arg20_l[8]);
    atomicAdd(&ind_arg6[9+map0idx*40],arg20_l[9]);
    atomicAdd(&ind_arg6[10+map0idx*40],arg20_l[10]);
    atomicAdd(&ind_arg6[11+map0idx*40],arg20_l[11]);
    atomicAdd(&ind_arg6[12+map0idx*40],arg20_l[12]);
    atomicAdd(&ind_arg6[13+map0idx*40],arg20_l[13]);
    atomicAdd(&ind_arg6[14+map0idx*40],arg20_l[14]);
    atomicAdd(&ind_arg6[15+map0idx*40],arg20_l[15]);
    atomicAdd(&ind_arg6[16+map0idx*40],arg20_l[16]);
    atomicAdd(&ind_arg6[17+map0idx*40],arg20_l[17]);
    atomicAdd(&ind_arg6[18+map0idx*40],arg20_l[18]);
    atomicAdd(&ind_arg6[19+map0idx*40],arg20_l[19]);
    atomicAdd(&ind_arg6[20+map0idx*40],arg20_l[20]);
    atomicAdd(&ind_arg6[21+map0idx*40],arg20_l[21]);
    atomicAdd(&ind_arg6[22+map0idx*40],arg20_l[22]);
    atomicAdd(&ind_arg6[23+map0idx*40],arg20_l[23]);
    atomicAdd(&ind_arg6[24+map0idx*40],arg20_l[24]);
    atomicAdd(&ind_arg6[25+map0idx*40],arg20_l[25]);
    atomicAdd(&ind_arg6[26+map0idx*40],arg20_l[26]);
    atomicAdd(&ind_arg6[27+map0idx*40],arg20_l[27]);
    atomicAdd(&ind_arg6[28+map0idx*40],arg20_l[28]);
    atomicAdd(&ind_arg6[29+map0idx*40],arg20_l[29]);
    atomicAdd(&ind_arg6[30+map0idx*40],arg20_l[30]);
    atomicAdd(&ind_arg6[31+map0idx*40],arg20_l[31]);
    atomicAdd(&ind_arg6[32+map0idx*40],arg20_l[32]);
    atomicAdd(&ind_arg6[33+map0idx*40],arg20_l[33]);
    atomicAdd(&ind_arg6[34+map0idx*40],arg20_l[34]);
    atomicAdd(&ind_arg6[35+map0idx*40],arg20_l[35]);
    atomicAdd(&ind_arg6[36+map0idx*40],arg20_l[36]);
    atomicAdd(&ind_arg6[37+map0idx*40],arg20_l[37]);
    atomicAdd(&ind_arg6[38+map0idx*40],arg20_l[38]);
    atomicAdd(&ind_arg6[39+map0idx*40],arg20_l[39]);
    atomicAdd(&ind_arg6[0+map1idx*40],arg21_l[0]);
    atomicAdd(&ind_arg6[1+map1idx*40],arg21_l[1]);
    atomicAdd(&ind_arg6[2+map1idx*40],arg21_l[2]);
    atomicAdd(&ind_arg6[3+map1idx*40],arg21_l[3]);
    atomicAdd(&ind_arg6[4+map1idx*40],arg21_l[4]);
    atomicAdd(&ind_arg6[5+map1idx*40],arg21_l[5]);
    atomicAdd(&ind_arg6[6+map1idx*40],arg21_l[6]);
    atomicAdd(&ind_arg6[7+map1idx*40],arg21_l[7]);
    atomicAdd(&ind_arg6[8+map1idx*40],arg21_l[8]);
    atomicAdd(&ind_arg6[9+map1idx*40],arg21_l[9]);
    atomicAdd(&ind_arg6[10+map1idx*40],arg21_l[10]);
    atomicAdd(&ind_arg6[11+map1idx*40],arg21_l[11]);
    atomicAdd(&ind_arg6[12+map1idx*40],arg21_l[12]);
    atomicAdd(&ind_arg6[13+map1idx*40],arg21_l[13]);
    atomicAdd(&ind_arg6[14+map1idx*40],arg21_l[14]);
    atomicAdd(&ind_arg6[15+map1idx*40],arg21_l[15]);
    atomicAdd(&ind_arg6[16+map1idx*40],arg21_l[16]);
    atomicAdd(&ind_arg6[17+map1idx*40],arg21_l[17]);
    atomicAdd(&ind_arg6[18+map1idx*40],arg21_l[18]);
    atomicAdd(&ind_arg6[19+map1idx*40],arg21_l[19]);
    atomicAdd(&ind_arg6[20+map1idx*40],arg21_l[20]);
    atomicAdd(&ind_arg6[21+map1idx*40],arg21_l[21]);
    atomicAdd(&ind_arg6[22+map1idx*40],arg21_l[22]);
    atomicAdd(&ind_arg6[23+map1idx*40],arg21_l[23]);
    atomicAdd(&ind_arg6[24+map1idx*40],arg21_l[24]);
    atomicAdd(&ind_arg6[25+map1idx*40],arg21_l[25]);
    atomicAdd(&ind_arg6[26+map1idx*40],arg21_l[26]);
    atomicAdd(&ind_arg6[27+map1idx*40],arg21_l[27]);
    atomicAdd(&ind_arg6[28+map1idx*40],arg21_l[28]);
    atomicAdd(&ind_arg6[29+map1idx*40],arg21_l[29]);
    atomicAdd(&ind_arg6[30+map1idx*40],arg21_l[30]);
    atomicAdd(&ind_arg6[31+map1idx*40],arg21_l[31]);
    atomicAdd(&ind_arg6[32+map1idx*40],arg21_l[32]);
    atomicAdd(&ind_arg6[33+map1idx*40],arg21_l[33]);
    atomicAdd(&ind_arg6[34+map1idx*40],arg21_l[34]);
    atomicAdd(&ind_arg6[35+map1idx*40],arg21_l[35]);
    atomicAdd(&ind_arg6[36+map1idx*40],arg21_l[36]);
    atomicAdd(&ind_arg6[37+map1idx*40],arg21_l[37]);
    atomicAdd(&ind_arg6[38+map1idx*40],arg21_l[38]);
    atomicAdd(&ind_arg6[39+map1idx*40],arg21_l[39]);
    atomicAdd(&ind_arg7[0+map0idx*40],arg22_l[0]);
    atomicAdd(&ind_arg7[1+map0idx*40],arg22_l[1]);
    atomicAdd(&ind_arg7[2+map0idx*40],arg22_l[2]);
    atomicAdd(&ind_arg7[3+map0idx*40],arg22_l[3]);
    atomicAdd(&ind_arg7[4+map0idx*40],arg22_l[4]);
    atomicAdd(&ind_arg7[5+map0idx*40],arg22_l[5]);
    atomicAdd(&ind_arg7[6+map0idx*40],arg22_l[6]);
    atomicAdd(&ind_arg7[7+map0idx*40],arg22_l[7]);
    atomicAdd(&ind_arg7[8+map0idx*40],arg22_l[8]);
    atomicAdd(&ind_arg7[9+map0idx*40],arg22_l[9]);
    atomicAdd(&ind_arg7[10+map0idx*40],arg22_l[10]);
    atomicAdd(&ind_arg7[11+map0idx*40],arg22_l[11]);
    atomicAdd(&ind_arg7[12+map0idx*40],arg22_l[12]);
    atomicAdd(&ind_arg7[13+map0idx*40],arg22_l[13]);
    atomicAdd(&ind_arg7[14+map0idx*40],arg22_l[14]);
    atomicAdd(&ind_arg7[15+map0idx*40],arg22_l[15]);
    atomicAdd(&ind_arg7[16+map0idx*40],arg22_l[16]);
    atomicAdd(&ind_arg7[17+map0idx*40],arg22_l[17]);
    atomicAdd(&ind_arg7[18+map0idx*40],arg22_l[18]);
    atomicAdd(&ind_arg7[19+map0idx*40],arg22_l[19]);
    atomicAdd(&ind_arg7[20+map0idx*40],arg22_l[20]);
    atomicAdd(&ind_arg7[21+map0idx*40],arg22_l[21]);
    atomicAdd(&ind_arg7[22+map0idx*40],arg22_l[22]);
    atomicAdd(&ind_arg7[23+map0idx*40],arg22_l[23]);
    atomicAdd(&ind_arg7[24+map0idx*40],arg22_l[24]);
    atomicAdd(&ind_arg7[25+map0idx*40],arg22_l[25]);
    atomicAdd(&ind_arg7[26+map0idx*40],arg22_l[26]);
    atomicAdd(&ind_arg7[27+map0idx*40],arg22_l[27]);
    atomicAdd(&ind_arg7[28+map0idx*40],arg22_l[28]);
    atomicAdd(&ind_arg7[29+map0idx*40],arg22_l[29]);
    atomicAdd(&ind_arg7[30+map0idx*40],arg22_l[30]);
    atomicAdd(&ind_arg7[31+map0idx*40],arg22_l[31]);
    atomicAdd(&ind_arg7[32+map0idx*40],arg22_l[32]);
    atomicAdd(&ind_arg7[33+map0idx*40],arg22_l[33]);
    atomicAdd(&ind_arg7[34+map0idx*40],arg22_l[34]);
    atomicAdd(&ind_arg7[35+map0idx*40],arg22_l[35]);
    atomicAdd(&ind_arg7[36+map0idx*40],arg22_l[36]);
    atomicAdd(&ind_arg7[37+map0idx*40],arg22_l[37]);
    atomicAdd(&ind_arg7[38+map0idx*40],arg22_l[38]);
    atomicAdd(&ind_arg7[39+map0idx*40],arg22_l[39]);
    atomicAdd(&ind_arg7[0+map1idx*40],arg23_l[0]);
    atomicAdd(&ind_arg7[1+map1idx*40],arg23_l[1]);
    atomicAdd(&ind_arg7[2+map1idx*40],arg23_l[2]);
    atomicAdd(&ind_arg7[3+map1idx*40],arg23_l[3]);
    atomicAdd(&ind_arg7[4+map1idx*40],arg23_l[4]);
    atomicAdd(&ind_arg7[5+map1idx*40],arg23_l[5]);
    atomicAdd(&ind_arg7[6+map1idx*40],arg23_l[6]);
    atomicAdd(&ind_arg7[7+map1idx*40],arg23_l[7]);
    atomicAdd(&ind_arg7[8+map1idx*40],arg23_l[8]);
    atomicAdd(&ind_arg7[9+map1idx*40],arg23_l[9]);
    atomicAdd(&ind_arg7[10+map1idx*40],arg23_l[10]);
    atomicAdd(&ind_arg7[11+map1idx*40],arg23_l[11]);
    atomicAdd(&ind_arg7[12+map1idx*40],arg23_l[12]);
    atomicAdd(&ind_arg7[13+map1idx*40],arg23_l[13]);
    atomicAdd(&ind_arg7[14+map1idx*40],arg23_l[14]);
    atomicAdd(&ind_arg7[15+map1idx*40],arg23_l[15]);
    atomicAdd(&ind_arg7[16+map1idx*40],arg23_l[16]);
    atomicAdd(&ind_arg7[17+map1idx*40],arg23_l[17]);
    atomicAdd(&ind_arg7[18+map1idx*40],arg23_l[18]);
    atomicAdd(&ind_arg7[19+map1idx*40],arg23_l[19]);
    atomicAdd(&ind_arg7[20+map1idx*40],arg23_l[20]);
    atomicAdd(&ind_arg7[21+map1idx*40],arg23_l[21]);
    atomicAdd(&ind_arg7[22+map1idx*40],arg23_l[22]);
    atomicAdd(&ind_arg7[23+map1idx*40],arg23_l[23]);
    atomicAdd(&ind_arg7[24+map1idx*40],arg23_l[24]);
    atomicAdd(&ind_arg7[25+map1idx*40],arg23_l[25]);
    atomicAdd(&ind_arg7[26+map1idx*40],arg23_l[26]);
    atomicAdd(&ind_arg7[27+map1idx*40],arg23_l[27]);
    atomicAdd(&ind_arg7[28+map1idx*40],arg23_l[28]);
    atomicAdd(&ind_arg7[29+map1idx*40],arg23_l[29]);
    atomicAdd(&ind_arg7[30+map1idx*40],arg23_l[30]);
    atomicAdd(&ind_arg7[31+map1idx*40],arg23_l[31]);
    atomicAdd(&ind_arg7[32+map1idx*40],arg23_l[32]);
    atomicAdd(&ind_arg7[33+map1idx*40],arg23_l[33]);
    atomicAdd(&ind_arg7[34+map1idx*40],arg23_l[34]);
    atomicAdd(&ind_arg7[35+map1idx*40],arg23_l[35]);
    atomicAdd(&ind_arg7[36+map1idx*40],arg23_l[36]);
    atomicAdd(&ind_arg7[37+map1idx*40],arg23_l[37]);
    atomicAdd(&ind_arg7[38+map1idx*40],arg23_l[38]);
    atomicAdd(&ind_arg7[39+map1idx*40],arg23_l[39]);
    atomicAdd(&ind_arg8[0+map0idx*40],arg24_l[0]);
    atomicAdd(&ind_arg8[1+map0idx*40],arg24_l[1]);
    atomicAdd(&ind_arg8[2+map0idx*40],arg24_l[2]);
    atomicAdd(&ind_arg8[3+map0idx*40],arg24_l[3]);
    atomicAdd(&ind_arg8[4+map0idx*40],arg24_l[4]);
    atomicAdd(&ind_arg8[5+map0idx*40],arg24_l[5]);
    atomicAdd(&ind_arg8[6+map0idx*40],arg24_l[6]);
    atomicAdd(&ind_arg8[7+map0idx*40],arg24_l[7]);
    atomicAdd(&ind_arg8[8+map0idx*40],arg24_l[8]);
    atomicAdd(&ind_arg8[9+map0idx*40],arg24_l[9]);
    atomicAdd(&ind_arg8[10+map0idx*40],arg24_l[10]);
    atomicAdd(&ind_arg8[11+map0idx*40],arg24_l[11]);
    atomicAdd(&ind_arg8[12+map0idx*40],arg24_l[12]);
    atomicAdd(&ind_arg8[13+map0idx*40],arg24_l[13]);
    atomicAdd(&ind_arg8[14+map0idx*40],arg24_l[14]);
    atomicAdd(&ind_arg8[15+map0idx*40],arg24_l[15]);
    atomicAdd(&ind_arg8[16+map0idx*40],arg24_l[16]);
    atomicAdd(&ind_arg8[17+map0idx*40],arg24_l[17]);
    atomicAdd(&ind_arg8[18+map0idx*40],arg24_l[18]);
    atomicAdd(&ind_arg8[19+map0idx*40],arg24_l[19]);
    atomicAdd(&ind_arg8[20+map0idx*40],arg24_l[20]);
    atomicAdd(&ind_arg8[21+map0idx*40],arg24_l[21]);
    atomicAdd(&ind_arg8[22+map0idx*40],arg24_l[22]);
    atomicAdd(&ind_arg8[23+map0idx*40],arg24_l[23]);
    atomicAdd(&ind_arg8[24+map0idx*40],arg24_l[24]);
    atomicAdd(&ind_arg8[25+map0idx*40],arg24_l[25]);
    atomicAdd(&ind_arg8[26+map0idx*40],arg24_l[26]);
    atomicAdd(&ind_arg8[27+map0idx*40],arg24_l[27]);
    atomicAdd(&ind_arg8[28+map0idx*40],arg24_l[28]);
    atomicAdd(&ind_arg8[29+map0idx*40],arg24_l[29]);
    atomicAdd(&ind_arg8[30+map0idx*40],arg24_l[30]);
    atomicAdd(&ind_arg8[31+map0idx*40],arg24_l[31]);
    atomicAdd(&ind_arg8[32+map0idx*40],arg24_l[32]);
    atomicAdd(&ind_arg8[33+map0idx*40],arg24_l[33]);
    atomicAdd(&ind_arg8[34+map0idx*40],arg24_l[34]);
    atomicAdd(&ind_arg8[35+map0idx*40],arg24_l[35]);
    atomicAdd(&ind_arg8[36+map0idx*40],arg24_l[36]);
    atomicAdd(&ind_arg8[37+map0idx*40],arg24_l[37]);
    atomicAdd(&ind_arg8[38+map0idx*40],arg24_l[38]);
    atomicAdd(&ind_arg8[39+map0idx*40],arg24_l[39]);
    atomicAdd(&ind_arg8[0+map1idx*40],arg25_l[0]);
    atomicAdd(&ind_arg8[1+map1idx*40],arg25_l[1]);
    atomicAdd(&ind_arg8[2+map1idx*40],arg25_l[2]);
    atomicAdd(&ind_arg8[3+map1idx*40],arg25_l[3]);
    atomicAdd(&ind_arg8[4+map1idx*40],arg25_l[4]);
    atomicAdd(&ind_arg8[5+map1idx*40],arg25_l[5]);
    atomicAdd(&ind_arg8[6+map1idx*40],arg25_l[6]);
    atomicAdd(&ind_arg8[7+map1idx*40],arg25_l[7]);
    atomicAdd(&ind_arg8[8+map1idx*40],arg25_l[8]);
    atomicAdd(&ind_arg8[9+map1idx*40],arg25_l[9]);
    atomicAdd(&ind_arg8[10+map1idx*40],arg25_l[10]);
    atomicAdd(&ind_arg8[11+map1idx*40],arg25_l[11]);
    atomicAdd(&ind_arg8[12+map1idx*40],arg25_l[12]);
    atomicAdd(&ind_arg8[13+map1idx*40],arg25_l[13]);
    atomicAdd(&ind_arg8[14+map1idx*40],arg25_l[14]);
    atomicAdd(&ind_arg8[15+map1idx*40],arg25_l[15]);
    atomicAdd(&ind_arg8[16+map1idx*40],arg25_l[16]);
    atomicAdd(&ind_arg8[17+map1idx*40],arg25_l[17]);
    atomicAdd(&ind_arg8[18+map1idx*40],arg25_l[18]);
    atomicAdd(&ind_arg8[19+map1idx*40],arg25_l[19]);
    atomicAdd(&ind_arg8[20+map1idx*40],arg25_l[20]);
    atomicAdd(&ind_arg8[21+map1idx*40],arg25_l[21]);
    atomicAdd(&ind_arg8[22+map1idx*40],arg25_l[22]);
    atomicAdd(&ind_arg8[23+map1idx*40],arg25_l[23]);
    atomicAdd(&ind_arg8[24+map1idx*40],arg25_l[24]);
    atomicAdd(&ind_arg8[25+map1idx*40],arg25_l[25]);
    atomicAdd(&ind_arg8[26+map1idx*40],arg25_l[26]);
    atomicAdd(&ind_arg8[27+map1idx*40],arg25_l[27]);
    atomicAdd(&ind_arg8[28+map1idx*40],arg25_l[28]);
    atomicAdd(&ind_arg8[29+map1idx*40],arg25_l[29]);
    atomicAdd(&ind_arg8[30+map1idx*40],arg25_l[30]);
    atomicAdd(&ind_arg8[31+map1idx*40],arg25_l[31]);
    atomicAdd(&ind_arg8[32+map1idx*40],arg25_l[32]);
    atomicAdd(&ind_arg8[33+map1idx*40],arg25_l[33]);
    atomicAdd(&ind_arg8[34+map1idx*40],arg25_l[34]);
    atomicAdd(&ind_arg8[35+map1idx*40],arg25_l[35]);
    atomicAdd(&ind_arg8[36+map1idx*40],arg25_l[36]);
    atomicAdd(&ind_arg8[37+map1idx*40],arg25_l[37]);
    atomicAdd(&ind_arg8[38+map1idx*40],arg25_l[38]);
    atomicAdd(&ind_arg8[39+map1idx*40],arg25_l[39]);
  }
}

//host stub function
void custom_kernel_pmf_3d_mult_faces(char const *name, op_set set,
  op_arg arg0,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg12,
  op_arg arg14,
  op_arg arg16,
  op_arg arg18,
  op_arg arg20,
  op_arg arg22,
  op_arg arg24){

  int nargs = 26;
  op_arg args[26];

  arg0.idx = 0;
  args[0] = arg0;
  for ( int v=1; v<2; v++ ){
    args[0 + v] = op_arg_dat(arg0.dat, v, arg0.map, 1, "int", OP_READ);
  }

  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 20, DG_FP_STR, OP_READ);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 20, DG_FP_STR, OP_READ);
  }

  arg14.idx = 0;
  args[14] = arg14;
  for ( int v=1; v<2; v++ ){
    args[14 + v] = op_arg_dat(arg14.dat, v, arg14.map, 20, DG_FP_STR, OP_READ);
  }

  arg16.idx = 0;
  args[16] = arg16;
  for ( int v=1; v<2; v++ ){
    args[16 + v] = op_arg_dat(arg16.dat, v, arg16.map, 20, DG_FP_STR, OP_READ);
  }

  arg18.idx = 0;
  args[18] = arg18;
  for ( int v=1; v<2; v++ ){
    args[18 + v] = op_arg_dat(arg18.dat, v, arg18.map, 40, DG_FP_STR, OP_INC);
  }

  arg20.idx = 0;
  args[20] = arg20;
  for ( int v=1; v<2; v++ ){
    args[20 + v] = op_arg_dat(arg20.dat, v, arg20.map, 40, DG_FP_STR, OP_INC);
  }

  arg22.idx = 0;
  args[22] = arg22;
  for ( int v=1; v<2; v++ ){
    args[22 + v] = op_arg_dat(arg22.dat, v, arg22.map, 40, DG_FP_STR, OP_INC);
  }

  arg24.idx = 0;
  args[24] = arg24;
  for ( int v=1; v<2; v++ ){
    args[24 + v] = op_arg_dat(arg24.dat, v, arg24.map, 40, DG_FP_STR, OP_INC);
  }

  int ninds   = 9;
  int inds[26] = {0,0,-1,-1,-1,-1,-1,-1,-1,-1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pmf_3d_mult_faces\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  if (set_size > 0) {

    //set CUDA execution parameters
    // #ifdef OP_BLOCK_SIZE_37
    //   int nthread = OP_BLOCK_SIZE_37;
    // #else
    //   int nthread = OP_block_size;
    // #endif
    const int nthread = 256;
    const int num_cells = (nthread / DG_NPF) + 2;

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        int nblocks = (end-start-1)/nthread+1;
        _op_cuda_pmf_3d_mult_faces<<<nblocks,nthread>>>(
        (int *)arg0.data_d,
        (DG_FP *)arg10.data_d,
        (DG_FP *)arg12.data_d,
        (DG_FP *)arg14.data_d,
        (DG_FP *)arg16.data_d,
        (DG_FP *)arg18.data_d,
        (DG_FP *)arg20.data_d,
        (DG_FP *)arg22.data_d,
        (DG_FP *)arg24.data_d,
        arg0.map_data_d,
        (int*)arg2.data_d,
        (int*)arg3.data_d,
        (int*)arg4.data_d,
        (DG_FP*)arg5.data_d,
        (DG_FP*)arg6.data_d,
        (DG_FP*)arg7.data_d,
        (DG_FP*)arg8.data_d,
        (DG_FP*)arg9.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
