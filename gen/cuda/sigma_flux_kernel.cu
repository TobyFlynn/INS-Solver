//
// auto-generated by op2.py
//

//user function
__device__ void sigma_flux_gpu( const int *edgeNum, const bool *rev, const double **sJ,
                       const double **nx, const double **ny, const double **s,
                       const double **vis, double **sigFx, double **sigFy) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exIndL = edgeL * 6;
  int exIndR = edgeR * 6;

  double kL = sqrt(vis[0][0]);
  double kR = sqrt(vis[1][0]);
  double lamdaL = kL;
  double lamdaR = kR;
  double wL, wR;

  if(lamdaL < 1e-12 && lamdaR < 1e-12) {
    wL = 0.5;
    wR = 0.5;
  } else {
    double lamdaAvg = (lamdaL + lamdaR) / 2.0;
    wL = lamdaL / (2.0 * lamdaAvg);
    wR = lamdaR / (2.0 * lamdaAvg);
  }

  for(int i = 0; i < 6; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + 6 - i - 1;
    } else {
      rInd = exIndR + i;
    }
    double flux = wL * s[0][lInd] + wR * s[1][rInd];
    flux *= kL;
    sigFx[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * nx[0][lInd] * flux;
    sigFy[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * ny[0][lInd] * flux;
  }

  for(int i = 0; i < 6; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 6 - i - 1;
    } else {
      lInd = exIndL + i;
    }
    double flux = wL * s[0][lInd] + wR * s[1][rInd];
    flux *= kR;
    sigFx[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * nx[1][rInd] * flux;
    sigFy[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * ny[1][rInd] * flux;
  }

}

// CUDA kernel function
__global__ void op_cuda_sigma_flux(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  double *__restrict ind_arg5,
  double *__restrict ind_arg6,
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
    double arg12_l[18];
    for ( int d=0; d<18; d++ ){
      arg12_l[d] = ZERO_double;
    }
    double arg13_l[18];
    for ( int d=0; d<18; d++ ){
      arg13_l[d] = ZERO_double;
    }
    double arg14_l[18];
    for ( int d=0; d<18; d++ ){
      arg14_l[d] = ZERO_double;
    }
    double arg15_l[18];
    for ( int d=0; d<18; d++ ){
      arg15_l[d] = ZERO_double;
    }
    int map2idx;
    int map3idx;
    map2idx = opDat2Map[n + set_size * 0];
    map3idx = opDat2Map[n + set_size * 1];
    const double* arg2_vec[] = {
       &ind_arg0[18 * map2idx],
       &ind_arg0[18 * map3idx]};
    const double* arg4_vec[] = {
       &ind_arg1[18 * map2idx],
       &ind_arg1[18 * map3idx]};
    const double* arg6_vec[] = {
       &ind_arg2[18 * map2idx],
       &ind_arg2[18 * map3idx]};
    const double* arg8_vec[] = {
       &ind_arg3[18 * map2idx],
       &ind_arg3[18 * map3idx]};
    const double* arg10_vec[] = {
       &ind_arg4[1 * map2idx],
       &ind_arg4[1 * map3idx]};
    double* arg12_vec[] = {
      arg12_l,
      arg13_l};
    double* arg14_vec[] = {
      arg14_l,
      arg15_l};

    //user-supplied kernel call
    sigma_flux_gpu(arg0+n*2,
               arg1+n*1,
               arg2_vec,
               arg4_vec,
               arg6_vec,
               arg8_vec,
               arg10_vec,
               arg12_vec,
               arg14_vec);
    atomicAdd(&ind_arg5[0+map2idx*18],arg12_l[0]);
    atomicAdd(&ind_arg5[1+map2idx*18],arg12_l[1]);
    atomicAdd(&ind_arg5[2+map2idx*18],arg12_l[2]);
    atomicAdd(&ind_arg5[3+map2idx*18],arg12_l[3]);
    atomicAdd(&ind_arg5[4+map2idx*18],arg12_l[4]);
    atomicAdd(&ind_arg5[5+map2idx*18],arg12_l[5]);
    atomicAdd(&ind_arg5[6+map2idx*18],arg12_l[6]);
    atomicAdd(&ind_arg5[7+map2idx*18],arg12_l[7]);
    atomicAdd(&ind_arg5[8+map2idx*18],arg12_l[8]);
    atomicAdd(&ind_arg5[9+map2idx*18],arg12_l[9]);
    atomicAdd(&ind_arg5[10+map2idx*18],arg12_l[10]);
    atomicAdd(&ind_arg5[11+map2idx*18],arg12_l[11]);
    atomicAdd(&ind_arg5[12+map2idx*18],arg12_l[12]);
    atomicAdd(&ind_arg5[13+map2idx*18],arg12_l[13]);
    atomicAdd(&ind_arg5[14+map2idx*18],arg12_l[14]);
    atomicAdd(&ind_arg5[15+map2idx*18],arg12_l[15]);
    atomicAdd(&ind_arg5[16+map2idx*18],arg12_l[16]);
    atomicAdd(&ind_arg5[17+map2idx*18],arg12_l[17]);
    atomicAdd(&ind_arg5[0+map3idx*18],arg13_l[0]);
    atomicAdd(&ind_arg5[1+map3idx*18],arg13_l[1]);
    atomicAdd(&ind_arg5[2+map3idx*18],arg13_l[2]);
    atomicAdd(&ind_arg5[3+map3idx*18],arg13_l[3]);
    atomicAdd(&ind_arg5[4+map3idx*18],arg13_l[4]);
    atomicAdd(&ind_arg5[5+map3idx*18],arg13_l[5]);
    atomicAdd(&ind_arg5[6+map3idx*18],arg13_l[6]);
    atomicAdd(&ind_arg5[7+map3idx*18],arg13_l[7]);
    atomicAdd(&ind_arg5[8+map3idx*18],arg13_l[8]);
    atomicAdd(&ind_arg5[9+map3idx*18],arg13_l[9]);
    atomicAdd(&ind_arg5[10+map3idx*18],arg13_l[10]);
    atomicAdd(&ind_arg5[11+map3idx*18],arg13_l[11]);
    atomicAdd(&ind_arg5[12+map3idx*18],arg13_l[12]);
    atomicAdd(&ind_arg5[13+map3idx*18],arg13_l[13]);
    atomicAdd(&ind_arg5[14+map3idx*18],arg13_l[14]);
    atomicAdd(&ind_arg5[15+map3idx*18],arg13_l[15]);
    atomicAdd(&ind_arg5[16+map3idx*18],arg13_l[16]);
    atomicAdd(&ind_arg5[17+map3idx*18],arg13_l[17]);
    atomicAdd(&ind_arg6[0+map2idx*18],arg14_l[0]);
    atomicAdd(&ind_arg6[1+map2idx*18],arg14_l[1]);
    atomicAdd(&ind_arg6[2+map2idx*18],arg14_l[2]);
    atomicAdd(&ind_arg6[3+map2idx*18],arg14_l[3]);
    atomicAdd(&ind_arg6[4+map2idx*18],arg14_l[4]);
    atomicAdd(&ind_arg6[5+map2idx*18],arg14_l[5]);
    atomicAdd(&ind_arg6[6+map2idx*18],arg14_l[6]);
    atomicAdd(&ind_arg6[7+map2idx*18],arg14_l[7]);
    atomicAdd(&ind_arg6[8+map2idx*18],arg14_l[8]);
    atomicAdd(&ind_arg6[9+map2idx*18],arg14_l[9]);
    atomicAdd(&ind_arg6[10+map2idx*18],arg14_l[10]);
    atomicAdd(&ind_arg6[11+map2idx*18],arg14_l[11]);
    atomicAdd(&ind_arg6[12+map2idx*18],arg14_l[12]);
    atomicAdd(&ind_arg6[13+map2idx*18],arg14_l[13]);
    atomicAdd(&ind_arg6[14+map2idx*18],arg14_l[14]);
    atomicAdd(&ind_arg6[15+map2idx*18],arg14_l[15]);
    atomicAdd(&ind_arg6[16+map2idx*18],arg14_l[16]);
    atomicAdd(&ind_arg6[17+map2idx*18],arg14_l[17]);
    atomicAdd(&ind_arg6[0+map3idx*18],arg15_l[0]);
    atomicAdd(&ind_arg6[1+map3idx*18],arg15_l[1]);
    atomicAdd(&ind_arg6[2+map3idx*18],arg15_l[2]);
    atomicAdd(&ind_arg6[3+map3idx*18],arg15_l[3]);
    atomicAdd(&ind_arg6[4+map3idx*18],arg15_l[4]);
    atomicAdd(&ind_arg6[5+map3idx*18],arg15_l[5]);
    atomicAdd(&ind_arg6[6+map3idx*18],arg15_l[6]);
    atomicAdd(&ind_arg6[7+map3idx*18],arg15_l[7]);
    atomicAdd(&ind_arg6[8+map3idx*18],arg15_l[8]);
    atomicAdd(&ind_arg6[9+map3idx*18],arg15_l[9]);
    atomicAdd(&ind_arg6[10+map3idx*18],arg15_l[10]);
    atomicAdd(&ind_arg6[11+map3idx*18],arg15_l[11]);
    atomicAdd(&ind_arg6[12+map3idx*18],arg15_l[12]);
    atomicAdd(&ind_arg6[13+map3idx*18],arg15_l[13]);
    atomicAdd(&ind_arg6[14+map3idx*18],arg15_l[14]);
    atomicAdd(&ind_arg6[15+map3idx*18],arg15_l[15]);
    atomicAdd(&ind_arg6[16+map3idx*18],arg15_l[16]);
    atomicAdd(&ind_arg6[17+map3idx*18],arg15_l[17]);
  }
}


//host stub function
void op_par_loop_sigma_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6,
  op_arg arg8,
  op_arg arg10,
  op_arg arg12,
  op_arg arg14){

  int nargs = 16;
  op_arg args[16];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 18, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 18, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 18, "double", OP_READ);
  }

  arg8.idx = 0;
  args[8] = arg8;
  for ( int v=1; v<2; v++ ){
    args[8 + v] = op_arg_dat(arg8.dat, v, arg8.map, 18, "double", OP_READ);
  }

  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 1, "double", OP_READ);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 18, "double", OP_INC);
  }

  arg14.idx = 0;
  args[14] = arg14;
  for ( int v=1; v<2; v++ ){
    args[14 + v] = op_arg_dat(arg14.dat, v, arg14.map, 18, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(26);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[26].name      = name;
  OP_kernels[26].count    += 1;


  int    ninds   = 7;
  int    inds[16] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: sigma_flux\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_26
      int nthread = OP_BLOCK_SIZE_26;
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
        op_cuda_sigma_flux<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg4.data_d,
        (double *)arg6.data_d,
        (double *)arg8.data_d,
        (double *)arg10.data_d,
        (double *)arg12.data_d,
        (double *)arg14.data_d,
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
  OP_kernels[26].time     += wall_t2 - wall_t1;
}
