//
// auto-generated by op2.py
//

//user function
__device__ void ls_flux_gpu( const int *edgeNum, const bool *rev, const double **sJ,
                    const double **nx, const double **ny, const double **s,
                    double **dsldx, double **dsrdx, double **dsldy,
                    double **dsrdy) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

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
      rInd = exIndR + 7 - i - 1;
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
      lInd = exIndL + 7 - i - 1;
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
  double *__restrict ind_arg4,
  double *__restrict ind_arg5,
  double *__restrict ind_arg6,
  double *__restrict ind_arg7,
  const int *__restrict opDat2Map,
  const int *__restrict arg0,
  const bool *__restrict arg1,
  int start,
  int end,
  int   set_size) {
  double arg10_l[21];
  double arg11_l[21];
  double arg12_l[21];
  double arg13_l[21];
  double arg14_l[21];
  double arg15_l[21];
  double arg16_l[21];
  double arg17_l[21];
  double *arg10_vec[2] = {
    arg10_l,
    arg11_l,
  };
  double *arg12_vec[2] = {
    arg12_l,
    arg13_l,
  };
  double *arg14_vec[2] = {
    arg14_l,
    arg15_l,
  };
  double *arg16_vec[2] = {
    arg16_l,
    arg17_l,
  };
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg10_l[21];
    for ( int d=0; d<21; d++ ){
      arg10_l[d] = ZERO_double;
    }
    double arg11_l[21];
    for ( int d=0; d<21; d++ ){
      arg11_l[d] = ZERO_double;
    }
    double arg12_l[21];
    for ( int d=0; d<21; d++ ){
      arg12_l[d] = ZERO_double;
    }
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
    int map2idx;
    int map3idx;
    map2idx = opDat2Map[n + set_size * 0];
    map3idx = opDat2Map[n + set_size * 1];
    const double* arg2_vec[] = {
       &ind_arg0[21 * map2idx],
       &ind_arg0[21 * map3idx]};
    const double* arg4_vec[] = {
       &ind_arg1[21 * map2idx],
       &ind_arg1[21 * map3idx]};
    const double* arg6_vec[] = {
       &ind_arg2[21 * map2idx],
       &ind_arg2[21 * map3idx]};
    const double* arg8_vec[] = {
       &ind_arg3[21 * map2idx],
       &ind_arg3[21 * map3idx]};
    double* arg10_vec[] = {
       &ind_arg4[21 * map2idx],
       &ind_arg4[21 * map3idx]};
    double* arg12_vec[] = {
       &ind_arg5[21 * map2idx],
       &ind_arg5[21 * map3idx]};
    double* arg14_vec[] = {
       &ind_arg6[21 * map2idx],
       &ind_arg6[21 * map3idx]};
    double* arg16_vec[] = {
       &ind_arg7[21 * map2idx],
       &ind_arg7[21 * map3idx]};

    //user-supplied kernel call
    ls_flux_gpu(arg0+n*2,
            arg1+n*1,
            arg2_vec,
            arg4_vec,
            arg6_vec,
            arg8_vec,
            arg10_vec,
            arg12_vec,
            arg14_vec,
            arg16_vec);
    atomicAdd(&ind_arg4[0+map2idx*21],arg10_l[0]);
    atomicAdd(&ind_arg4[1+map2idx*21],arg10_l[1]);
    atomicAdd(&ind_arg4[2+map2idx*21],arg10_l[2]);
    atomicAdd(&ind_arg4[3+map2idx*21],arg10_l[3]);
    atomicAdd(&ind_arg4[4+map2idx*21],arg10_l[4]);
    atomicAdd(&ind_arg4[5+map2idx*21],arg10_l[5]);
    atomicAdd(&ind_arg4[6+map2idx*21],arg10_l[6]);
    atomicAdd(&ind_arg4[7+map2idx*21],arg10_l[7]);
    atomicAdd(&ind_arg4[8+map2idx*21],arg10_l[8]);
    atomicAdd(&ind_arg4[9+map2idx*21],arg10_l[9]);
    atomicAdd(&ind_arg4[10+map2idx*21],arg10_l[10]);
    atomicAdd(&ind_arg4[11+map2idx*21],arg10_l[11]);
    atomicAdd(&ind_arg4[12+map2idx*21],arg10_l[12]);
    atomicAdd(&ind_arg4[13+map2idx*21],arg10_l[13]);
    atomicAdd(&ind_arg4[14+map2idx*21],arg10_l[14]);
    atomicAdd(&ind_arg4[15+map2idx*21],arg10_l[15]);
    atomicAdd(&ind_arg4[16+map2idx*21],arg10_l[16]);
    atomicAdd(&ind_arg4[17+map2idx*21],arg10_l[17]);
    atomicAdd(&ind_arg4[18+map2idx*21],arg10_l[18]);
    atomicAdd(&ind_arg4[19+map2idx*21],arg10_l[19]);
    atomicAdd(&ind_arg4[20+map2idx*21],arg10_l[20]);
    atomicAdd(&ind_arg4[0+map3idx*21],arg11_l[0]);
    atomicAdd(&ind_arg4[1+map3idx*21],arg11_l[1]);
    atomicAdd(&ind_arg4[2+map3idx*21],arg11_l[2]);
    atomicAdd(&ind_arg4[3+map3idx*21],arg11_l[3]);
    atomicAdd(&ind_arg4[4+map3idx*21],arg11_l[4]);
    atomicAdd(&ind_arg4[5+map3idx*21],arg11_l[5]);
    atomicAdd(&ind_arg4[6+map3idx*21],arg11_l[6]);
    atomicAdd(&ind_arg4[7+map3idx*21],arg11_l[7]);
    atomicAdd(&ind_arg4[8+map3idx*21],arg11_l[8]);
    atomicAdd(&ind_arg4[9+map3idx*21],arg11_l[9]);
    atomicAdd(&ind_arg4[10+map3idx*21],arg11_l[10]);
    atomicAdd(&ind_arg4[11+map3idx*21],arg11_l[11]);
    atomicAdd(&ind_arg4[12+map3idx*21],arg11_l[12]);
    atomicAdd(&ind_arg4[13+map3idx*21],arg11_l[13]);
    atomicAdd(&ind_arg4[14+map3idx*21],arg11_l[14]);
    atomicAdd(&ind_arg4[15+map3idx*21],arg11_l[15]);
    atomicAdd(&ind_arg4[16+map3idx*21],arg11_l[16]);
    atomicAdd(&ind_arg4[17+map3idx*21],arg11_l[17]);
    atomicAdd(&ind_arg4[18+map3idx*21],arg11_l[18]);
    atomicAdd(&ind_arg4[19+map3idx*21],arg11_l[19]);
    atomicAdd(&ind_arg4[20+map3idx*21],arg11_l[20]);
    atomicAdd(&ind_arg5[0+map2idx*21],arg12_l[0]);
    atomicAdd(&ind_arg5[1+map2idx*21],arg12_l[1]);
    atomicAdd(&ind_arg5[2+map2idx*21],arg12_l[2]);
    atomicAdd(&ind_arg5[3+map2idx*21],arg12_l[3]);
    atomicAdd(&ind_arg5[4+map2idx*21],arg12_l[4]);
    atomicAdd(&ind_arg5[5+map2idx*21],arg12_l[5]);
    atomicAdd(&ind_arg5[6+map2idx*21],arg12_l[6]);
    atomicAdd(&ind_arg5[7+map2idx*21],arg12_l[7]);
    atomicAdd(&ind_arg5[8+map2idx*21],arg12_l[8]);
    atomicAdd(&ind_arg5[9+map2idx*21],arg12_l[9]);
    atomicAdd(&ind_arg5[10+map2idx*21],arg12_l[10]);
    atomicAdd(&ind_arg5[11+map2idx*21],arg12_l[11]);
    atomicAdd(&ind_arg5[12+map2idx*21],arg12_l[12]);
    atomicAdd(&ind_arg5[13+map2idx*21],arg12_l[13]);
    atomicAdd(&ind_arg5[14+map2idx*21],arg12_l[14]);
    atomicAdd(&ind_arg5[15+map2idx*21],arg12_l[15]);
    atomicAdd(&ind_arg5[16+map2idx*21],arg12_l[16]);
    atomicAdd(&ind_arg5[17+map2idx*21],arg12_l[17]);
    atomicAdd(&ind_arg5[18+map2idx*21],arg12_l[18]);
    atomicAdd(&ind_arg5[19+map2idx*21],arg12_l[19]);
    atomicAdd(&ind_arg5[20+map2idx*21],arg12_l[20]);
    atomicAdd(&ind_arg5[0+map3idx*21],arg13_l[0]);
    atomicAdd(&ind_arg5[1+map3idx*21],arg13_l[1]);
    atomicAdd(&ind_arg5[2+map3idx*21],arg13_l[2]);
    atomicAdd(&ind_arg5[3+map3idx*21],arg13_l[3]);
    atomicAdd(&ind_arg5[4+map3idx*21],arg13_l[4]);
    atomicAdd(&ind_arg5[5+map3idx*21],arg13_l[5]);
    atomicAdd(&ind_arg5[6+map3idx*21],arg13_l[6]);
    atomicAdd(&ind_arg5[7+map3idx*21],arg13_l[7]);
    atomicAdd(&ind_arg5[8+map3idx*21],arg13_l[8]);
    atomicAdd(&ind_arg5[9+map3idx*21],arg13_l[9]);
    atomicAdd(&ind_arg5[10+map3idx*21],arg13_l[10]);
    atomicAdd(&ind_arg5[11+map3idx*21],arg13_l[11]);
    atomicAdd(&ind_arg5[12+map3idx*21],arg13_l[12]);
    atomicAdd(&ind_arg5[13+map3idx*21],arg13_l[13]);
    atomicAdd(&ind_arg5[14+map3idx*21],arg13_l[14]);
    atomicAdd(&ind_arg5[15+map3idx*21],arg13_l[15]);
    atomicAdd(&ind_arg5[16+map3idx*21],arg13_l[16]);
    atomicAdd(&ind_arg5[17+map3idx*21],arg13_l[17]);
    atomicAdd(&ind_arg5[18+map3idx*21],arg13_l[18]);
    atomicAdd(&ind_arg5[19+map3idx*21],arg13_l[19]);
    atomicAdd(&ind_arg5[20+map3idx*21],arg13_l[20]);
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
    atomicAdd(&ind_arg6[0+map3idx*21],arg15_l[0]);
    atomicAdd(&ind_arg6[1+map3idx*21],arg15_l[1]);
    atomicAdd(&ind_arg6[2+map3idx*21],arg15_l[2]);
    atomicAdd(&ind_arg6[3+map3idx*21],arg15_l[3]);
    atomicAdd(&ind_arg6[4+map3idx*21],arg15_l[4]);
    atomicAdd(&ind_arg6[5+map3idx*21],arg15_l[5]);
    atomicAdd(&ind_arg6[6+map3idx*21],arg15_l[6]);
    atomicAdd(&ind_arg6[7+map3idx*21],arg15_l[7]);
    atomicAdd(&ind_arg6[8+map3idx*21],arg15_l[8]);
    atomicAdd(&ind_arg6[9+map3idx*21],arg15_l[9]);
    atomicAdd(&ind_arg6[10+map3idx*21],arg15_l[10]);
    atomicAdd(&ind_arg6[11+map3idx*21],arg15_l[11]);
    atomicAdd(&ind_arg6[12+map3idx*21],arg15_l[12]);
    atomicAdd(&ind_arg6[13+map3idx*21],arg15_l[13]);
    atomicAdd(&ind_arg6[14+map3idx*21],arg15_l[14]);
    atomicAdd(&ind_arg6[15+map3idx*21],arg15_l[15]);
    atomicAdd(&ind_arg6[16+map3idx*21],arg15_l[16]);
    atomicAdd(&ind_arg6[17+map3idx*21],arg15_l[17]);
    atomicAdd(&ind_arg6[18+map3idx*21],arg15_l[18]);
    atomicAdd(&ind_arg6[19+map3idx*21],arg15_l[19]);
    atomicAdd(&ind_arg6[20+map3idx*21],arg15_l[20]);
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
    atomicAdd(&ind_arg7[0+map3idx*21],arg17_l[0]);
    atomicAdd(&ind_arg7[1+map3idx*21],arg17_l[1]);
    atomicAdd(&ind_arg7[2+map3idx*21],arg17_l[2]);
    atomicAdd(&ind_arg7[3+map3idx*21],arg17_l[3]);
    atomicAdd(&ind_arg7[4+map3idx*21],arg17_l[4]);
    atomicAdd(&ind_arg7[5+map3idx*21],arg17_l[5]);
    atomicAdd(&ind_arg7[6+map3idx*21],arg17_l[6]);
    atomicAdd(&ind_arg7[7+map3idx*21],arg17_l[7]);
    atomicAdd(&ind_arg7[8+map3idx*21],arg17_l[8]);
    atomicAdd(&ind_arg7[9+map3idx*21],arg17_l[9]);
    atomicAdd(&ind_arg7[10+map3idx*21],arg17_l[10]);
    atomicAdd(&ind_arg7[11+map3idx*21],arg17_l[11]);
    atomicAdd(&ind_arg7[12+map3idx*21],arg17_l[12]);
    atomicAdd(&ind_arg7[13+map3idx*21],arg17_l[13]);
    atomicAdd(&ind_arg7[14+map3idx*21],arg17_l[14]);
    atomicAdd(&ind_arg7[15+map3idx*21],arg17_l[15]);
    atomicAdd(&ind_arg7[16+map3idx*21],arg17_l[16]);
    atomicAdd(&ind_arg7[17+map3idx*21],arg17_l[17]);
    atomicAdd(&ind_arg7[18+map3idx*21],arg17_l[18]);
    atomicAdd(&ind_arg7[19+map3idx*21],arg17_l[19]);
    atomicAdd(&ind_arg7[20+map3idx*21],arg17_l[20]);
  }
}


//host stub function
void op_par_loop_ls_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6,
  op_arg arg8,
  op_arg arg10,
  op_arg arg12,
  op_arg arg14,
  op_arg arg16){

  int nargs = 18;
  op_arg args[18];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 21, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 21, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 21, "double", OP_READ);
  }

  arg8.idx = 0;
  args[8] = arg8;
  for ( int v=1; v<2; v++ ){
    args[8 + v] = op_arg_dat(arg8.dat, v, arg8.map, 21, "double", OP_READ);
  }

  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 21, "double", OP_INC);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 21, "double", OP_INC);
  }

  arg14.idx = 0;
  args[14] = arg14;
  for ( int v=1; v<2; v++ ){
    args[14 + v] = op_arg_dat(arg14.dat, v, arg14.map, 21, "double", OP_INC);
  }

  arg16.idx = 0;
  args[16] = arg16;
  for ( int v=1; v<2; v++ ){
    args[16 + v] = op_arg_dat(arg16.dat, v, arg16.map, 21, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(62);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[62].name      = name;
  OP_kernels[62].count    += 1;


  int    ninds   = 8;
  int    inds[18] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_flux\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_62
      int nthread = OP_BLOCK_SIZE_62;
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
        (double *)arg2.data_d,
        (double *)arg4.data_d,
        (double *)arg6.data_d,
        (double *)arg8.data_d,
        (double *)arg10.data_d,
        (double *)arg12.data_d,
        (double *)arg14.data_d,
        (double *)arg16.data_d,
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
  OP_kernels[62].time     += wall_t2 - wall_t1;
}
