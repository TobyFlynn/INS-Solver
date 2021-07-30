//
// auto-generated by op2.py
//

//user function
__device__ void diff_flux_gpu( const int *edgeNum, const bool *rev, const double **sJ,
                      const double **nx, const double **ny, const double **sigX,
                      const double **sigY, double **flux) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exIndL = edgeL * 3;
  int exIndR = edgeR * 3;

  for(int i = 0; i < 3; i++) {
    int rInd;
    int lInd = exIndL + i;
    if(reverse) {
      rInd = exIndR + 3 - i - 1;
    } else {
      rInd = exIndR + i;
    }

    double sigFX = (sigX[0][lInd] + sigX[1][rInd]) / 2.0;
    double sigFY = (sigY[0][lInd] + sigY[1][rInd]) / 2.0;

    flux[0][lInd] += gaussW_g_cuda[i] * sJ[0][lInd] * (nx[0][lInd] * sigFX + ny[0][lInd] * sigFY);
  }

  for(int i = 0; i < 3; i++) {
    int lInd;
    int rInd = exIndR + i;
    if(reverse) {
      lInd = exIndL + 3 - i - 1;
    } else {
      lInd = exIndL + i;
    }

    double sigFX = (sigX[0][lInd] + sigX[1][rInd]) / 2.0;
    double sigFY = (sigY[0][lInd] + sigY[1][rInd]) / 2.0;

    flux[1][rInd] += gaussW_g_cuda[i] * sJ[1][rInd] * (nx[1][rInd] * sigFX + ny[1][rInd] * sigFY);
  }

}

// CUDA kernel function
__global__ void op_cuda_diff_flux(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  double *__restrict ind_arg5,
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
    double arg12_l[9];
    for ( int d=0; d<9; d++ ){
      arg12_l[d] = ZERO_double;
    }
    double arg13_l[9];
    for ( int d=0; d<9; d++ ){
      arg13_l[d] = ZERO_double;
    }
    int map2idx;
    int map3idx;
    map2idx = opDat2Map[n + set_size * 0];
    map3idx = opDat2Map[n + set_size * 1];
    const double* arg2_vec[] = {
       &ind_arg0[9 * map2idx],
       &ind_arg0[9 * map3idx]};
    const double* arg4_vec[] = {
       &ind_arg1[9 * map2idx],
       &ind_arg1[9 * map3idx]};
    const double* arg6_vec[] = {
       &ind_arg2[9 * map2idx],
       &ind_arg2[9 * map3idx]};
    const double* arg8_vec[] = {
       &ind_arg3[9 * map2idx],
       &ind_arg3[9 * map3idx]};
    const double* arg10_vec[] = {
       &ind_arg4[9 * map2idx],
       &ind_arg4[9 * map3idx]};
    double* arg12_vec[] = {
      arg12_l,
      arg13_l};

    //user-supplied kernel call
    diff_flux_gpu(arg0+n*2,
              arg1+n*1,
              arg2_vec,
              arg4_vec,
              arg6_vec,
              arg8_vec,
              arg10_vec,
              arg12_vec);
    atomicAdd(&ind_arg5[0+map2idx*9],arg12_l[0]);
    atomicAdd(&ind_arg5[1+map2idx*9],arg12_l[1]);
    atomicAdd(&ind_arg5[2+map2idx*9],arg12_l[2]);
    atomicAdd(&ind_arg5[3+map2idx*9],arg12_l[3]);
    atomicAdd(&ind_arg5[4+map2idx*9],arg12_l[4]);
    atomicAdd(&ind_arg5[5+map2idx*9],arg12_l[5]);
    atomicAdd(&ind_arg5[6+map2idx*9],arg12_l[6]);
    atomicAdd(&ind_arg5[7+map2idx*9],arg12_l[7]);
    atomicAdd(&ind_arg5[8+map2idx*9],arg12_l[8]);
    atomicAdd(&ind_arg5[0+map3idx*9],arg13_l[0]);
    atomicAdd(&ind_arg5[1+map3idx*9],arg13_l[1]);
    atomicAdd(&ind_arg5[2+map3idx*9],arg13_l[2]);
    atomicAdd(&ind_arg5[3+map3idx*9],arg13_l[3]);
    atomicAdd(&ind_arg5[4+map3idx*9],arg13_l[4]);
    atomicAdd(&ind_arg5[5+map3idx*9],arg13_l[5]);
    atomicAdd(&ind_arg5[6+map3idx*9],arg13_l[6]);
    atomicAdd(&ind_arg5[7+map3idx*9],arg13_l[7]);
    atomicAdd(&ind_arg5[8+map3idx*9],arg13_l[8]);
  }
}


//host stub function
void op_par_loop_diff_flux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6,
  op_arg arg8,
  op_arg arg10,
  op_arg arg12){

  int nargs = 14;
  op_arg args[14];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 9, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 9, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 9, "double", OP_READ);
  }

  arg8.idx = 0;
  args[8] = arg8;
  for ( int v=1; v<2; v++ ){
    args[8 + v] = op_arg_dat(arg8.dat, v, arg8.map, 9, "double", OP_READ);
  }

  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 9, "double", OP_READ);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 9, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(61);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[61].name      = name;
  OP_kernels[61].count    += 1;


  int    ninds   = 6;
  int    inds[14] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: diff_flux\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_61
      int nthread = OP_BLOCK_SIZE_61;
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
        op_cuda_diff_flux<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg4.data_d,
        (double *)arg6.data_d,
        (double *)arg8.data_d,
        (double *)arg10.data_d,
        (double *)arg12.data_d,
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
  OP_kernels[61].time     += wall_t2 - wall_t1;
}
