//
// auto-generated by op2.py
//

//user function
__device__ void set_tau_gpu( const int *edgeNum, const double **x, const double **y,
                    const double **J, const double **sJ, double **tau) {

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
  if(edgeL == 1) exIndL = 5;
  else if(edgeL == 2) exIndL = 2 * 5;

  int exIndR = 0;
  if(edgeR == 1) exIndR = 5;
  else if(edgeR == 2) exIndR = 2 * 5;

  int *fmaskR;

  if(edgeR == 0) {
    fmaskR = FMASK_cuda;
  } else if(edgeR == 1) {
    fmaskR = &FMASK_cuda[5];
  } else {
    fmaskR = &FMASK_cuda[2 * 5];
  }

  int *fmaskL;

  if(edgeL == 0) {
    fmaskL = FMASK_cuda;
  } else if(edgeL == 1) {
    fmaskL = &FMASK_cuda[5];
  } else {
    fmaskL = &FMASK_cuda[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int rIndF, lIndF, rInd, lInd;
    if(reverse) {
      rIndF = fmaskR[5 - i - 1];
      rInd = exIndR + 5 - i - 1;
    } else {
      rIndF = fmaskR[i];
      rInd = exIndR + i;
    }
    lIndF = fmaskL[i];
    lInd = exIndL + i;

    double lH = 2.0 * J[0][lIndF] / sJ[0][lInd];
    double rH = 2.0 * J[1][rIndF] / sJ[1][rInd];
    if(lH < rH) {
      tau[0][lInd] += 15.0 / lH;
      tau[1][rInd] += 15.0 / lH;
    } else {
      tau[0][lInd] += 15.0 / rH;
      tau[1][rInd] += 15.0 / rH;
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_set_tau(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  double *__restrict ind_arg4,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  int start,
  int end,
  int   set_size) {
  double arg9_l[15];
  double arg10_l[15];
  double *arg9_vec[2] = {
    arg9_l,
    arg10_l,
  };
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg9_l[15];
    for ( int d=0; d<15; d++ ){
      arg9_l[d] = ZERO_double;
    }
    double arg10_l[15];
    for ( int d=0; d<15; d++ ){
      arg10_l[d] = ZERO_double;
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
       &ind_arg2[15 * map1idx],
       &ind_arg2[15 * map2idx]};
    const double* arg7_vec[] = {
       &ind_arg3[15 * map1idx],
       &ind_arg3[15 * map2idx]};
    double* arg9_vec[] = {
       &ind_arg4[15 * map1idx],
       &ind_arg4[15 * map2idx]};

    //user-supplied kernel call
    set_tau_gpu(arg0+n*2,
            arg1_vec,
            arg3_vec,
            arg5_vec,
            arg7_vec,
            arg9_vec);
    atomicAdd(&ind_arg4[0+map1idx*15],arg9_l[0]);
    atomicAdd(&ind_arg4[1+map1idx*15],arg9_l[1]);
    atomicAdd(&ind_arg4[2+map1idx*15],arg9_l[2]);
    atomicAdd(&ind_arg4[3+map1idx*15],arg9_l[3]);
    atomicAdd(&ind_arg4[4+map1idx*15],arg9_l[4]);
    atomicAdd(&ind_arg4[5+map1idx*15],arg9_l[5]);
    atomicAdd(&ind_arg4[6+map1idx*15],arg9_l[6]);
    atomicAdd(&ind_arg4[7+map1idx*15],arg9_l[7]);
    atomicAdd(&ind_arg4[8+map1idx*15],arg9_l[8]);
    atomicAdd(&ind_arg4[9+map1idx*15],arg9_l[9]);
    atomicAdd(&ind_arg4[10+map1idx*15],arg9_l[10]);
    atomicAdd(&ind_arg4[11+map1idx*15],arg9_l[11]);
    atomicAdd(&ind_arg4[12+map1idx*15],arg9_l[12]);
    atomicAdd(&ind_arg4[13+map1idx*15],arg9_l[13]);
    atomicAdd(&ind_arg4[14+map1idx*15],arg9_l[14]);
    atomicAdd(&ind_arg4[0+map2idx*15],arg10_l[0]);
    atomicAdd(&ind_arg4[1+map2idx*15],arg10_l[1]);
    atomicAdd(&ind_arg4[2+map2idx*15],arg10_l[2]);
    atomicAdd(&ind_arg4[3+map2idx*15],arg10_l[3]);
    atomicAdd(&ind_arg4[4+map2idx*15],arg10_l[4]);
    atomicAdd(&ind_arg4[5+map2idx*15],arg10_l[5]);
    atomicAdd(&ind_arg4[6+map2idx*15],arg10_l[6]);
    atomicAdd(&ind_arg4[7+map2idx*15],arg10_l[7]);
    atomicAdd(&ind_arg4[8+map2idx*15],arg10_l[8]);
    atomicAdd(&ind_arg4[9+map2idx*15],arg10_l[9]);
    atomicAdd(&ind_arg4[10+map2idx*15],arg10_l[10]);
    atomicAdd(&ind_arg4[11+map2idx*15],arg10_l[11]);
    atomicAdd(&ind_arg4[12+map2idx*15],arg10_l[12]);
    atomicAdd(&ind_arg4[13+map2idx*15],arg10_l[13]);
    atomicAdd(&ind_arg4[14+map2idx*15],arg10_l[14]);
  }
}


//host stub function
void op_par_loop_set_tau(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3,
  op_arg arg5,
  op_arg arg7,
  op_arg arg9){

  int nargs = 11;
  op_arg args[11];

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
    args[5 + v] = op_arg_dat(arg5.dat, v, arg5.map, 15, "double", OP_READ);
  }

  arg7.idx = 0;
  args[7] = arg7;
  for ( int v=1; v<2; v++ ){
    args[7 + v] = op_arg_dat(arg7.dat, v, arg7.map, 15, "double", OP_READ);
  }

  arg9.idx = 0;
  args[9] = arg9;
  for ( int v=1; v<2; v++ ){
    args[9 + v] = op_arg_dat(arg9.dat, v, arg9.map, 15, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(10);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[10].name      = name;
  OP_kernels[10].count    += 1;


  int    ninds   = 5;
  int    inds[11] = {-1,0,0,1,1,2,2,3,3,4,4};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: set_tau\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_10
      int nthread = OP_BLOCK_SIZE_10;
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
        op_cuda_set_tau<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg3.data_d,
        (double *)arg5.data_d,
        (double *)arg7.data_d,
        (double *)arg9.data_d,
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
  OP_kernels[10].time     += wall_t2 - wall_t1;
}
