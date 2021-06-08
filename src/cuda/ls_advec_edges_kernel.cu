//
// auto-generated by op2.py
//

//user function
__device__ void ls_advec_edges_gpu( const int *edgeNum, const double **x,
                        const double **y, const double **q, double **exQ) {

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

  int exInd = 0;
  if(edgeL == 1) exInd = 5;
  else if(edgeL == 2) exInd = 2 * 5;

  int *fmask;

  if(edgeR == 0) {
    fmask = FMASK_cuda;
  } else if(edgeR == 1) {
    fmask = &FMASK_cuda[5];
  } else {
    fmask = &FMASK_cuda[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int rInd;
    if(reverse) {
      rInd = fmask[5 - i - 1];
    } else {
      rInd = fmask[i];
    }
    exQ[0][exInd + i] += q[1][rInd];
  }

  exInd = 0;
  if(edgeR == 1) exInd = 5;
  else if(edgeR == 2) exInd = 2 * 5;

  if(edgeL == 0) {
    fmask = FMASK_cuda;
  } else if(edgeL == 1) {
    fmask = &FMASK_cuda[5];
  } else {
    fmask = &FMASK_cuda[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    int lInd;
    if(reverse) {
      lInd = fmask[5 - i - 1];
    } else {
      lInd = fmask[i];
    }
    exQ[1][exInd + i] += q[0][lInd];
  }

}

// CUDA kernel function
__global__ void op_cuda_ls_advec_edges(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  double *__restrict ind_arg3,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  int start,
  int end,
  int   set_size) {
  double arg7_l[15];
  double arg8_l[15];
  double *arg7_vec[2] = {
    arg7_l,
    arg8_l,
  };
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg7_l[15];
    for ( int d=0; d<15; d++ ){
      arg7_l[d] = ZERO_double;
    }
    double arg8_l[15];
    for ( int d=0; d<15; d++ ){
      arg8_l[d] = ZERO_double;
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
    double* arg7_vec[] = {
       &ind_arg3[15 * map1idx],
       &ind_arg3[15 * map2idx]};

    //user-supplied kernel call
    ls_advec_edges_gpu(arg0+n*2,
                   arg1_vec,
                   arg3_vec,
                   arg5_vec,
                   arg7_vec);
    atomicAdd(&ind_arg3[0+map1idx*15],arg7_l[0]);
    atomicAdd(&ind_arg3[1+map1idx*15],arg7_l[1]);
    atomicAdd(&ind_arg3[2+map1idx*15],arg7_l[2]);
    atomicAdd(&ind_arg3[3+map1idx*15],arg7_l[3]);
    atomicAdd(&ind_arg3[4+map1idx*15],arg7_l[4]);
    atomicAdd(&ind_arg3[5+map1idx*15],arg7_l[5]);
    atomicAdd(&ind_arg3[6+map1idx*15],arg7_l[6]);
    atomicAdd(&ind_arg3[7+map1idx*15],arg7_l[7]);
    atomicAdd(&ind_arg3[8+map1idx*15],arg7_l[8]);
    atomicAdd(&ind_arg3[9+map1idx*15],arg7_l[9]);
    atomicAdd(&ind_arg3[10+map1idx*15],arg7_l[10]);
    atomicAdd(&ind_arg3[11+map1idx*15],arg7_l[11]);
    atomicAdd(&ind_arg3[12+map1idx*15],arg7_l[12]);
    atomicAdd(&ind_arg3[13+map1idx*15],arg7_l[13]);
    atomicAdd(&ind_arg3[14+map1idx*15],arg7_l[14]);
    atomicAdd(&ind_arg3[0+map2idx*15],arg8_l[0]);
    atomicAdd(&ind_arg3[1+map2idx*15],arg8_l[1]);
    atomicAdd(&ind_arg3[2+map2idx*15],arg8_l[2]);
    atomicAdd(&ind_arg3[3+map2idx*15],arg8_l[3]);
    atomicAdd(&ind_arg3[4+map2idx*15],arg8_l[4]);
    atomicAdd(&ind_arg3[5+map2idx*15],arg8_l[5]);
    atomicAdd(&ind_arg3[6+map2idx*15],arg8_l[6]);
    atomicAdd(&ind_arg3[7+map2idx*15],arg8_l[7]);
    atomicAdd(&ind_arg3[8+map2idx*15],arg8_l[8]);
    atomicAdd(&ind_arg3[9+map2idx*15],arg8_l[9]);
    atomicAdd(&ind_arg3[10+map2idx*15],arg8_l[10]);
    atomicAdd(&ind_arg3[11+map2idx*15],arg8_l[11]);
    atomicAdd(&ind_arg3[12+map2idx*15],arg8_l[12]);
    atomicAdd(&ind_arg3[13+map2idx*15],arg8_l[13]);
    atomicAdd(&ind_arg3[14+map2idx*15],arg8_l[14]);
  }
}


//host stub function
void op_par_loop_ls_advec_edges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3,
  op_arg arg5,
  op_arg arg7){

  int nargs = 9;
  op_arg args[9];

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
    args[7 + v] = op_arg_dat(arg7.dat, v, arg7.map, 15, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(55);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[55].name      = name;
  OP_kernels[55].count    += 1;


  int    ninds   = 4;
  int    inds[9] = {-1,0,0,1,1,2,2,3,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_advec_edges\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_55
      int nthread = OP_BLOCK_SIZE_55;
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
        op_cuda_ls_advec_edges<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg3.data_d,
        (double *)arg5.data_d,
        (double *)arg7.data_d,
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
  OP_kernels[55].time     += wall_t2 - wall_t1;
}
