//
// auto-generated by op2.py
//

//user function
__device__ void ls_advec_edges_gpu( const int *edgeNum, const bool *rev,
                           const double **q, double **exQ) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  int exInd = edgeL * 4;
  int *fmask = &FMASK_cuda[edgeR * 4];

  for(int i = 0; i < 4; i++) {
    int rInd;
    if(reverse) {
      rInd = fmask[4 - i - 1];
    } else {
      rInd = fmask[i];
    }
    exQ[0][exInd + i] += q[1][rInd];
  }

  exInd = edgeR * 4;
  fmask = &FMASK_cuda[edgeL * 4];

  for(int i = 0; i < 4; i++) {
    int lInd;
    if(reverse) {
      lInd = fmask[4 - i - 1];
    } else {
      lInd = fmask[i];
    }
    exQ[1][exInd + i] += q[0][lInd];
  }

}

// CUDA kernel function
__global__ void op_cuda_ls_advec_edges(
  const double *__restrict ind_arg0,
  double *__restrict ind_arg1,
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
    double arg4_l[12];
    for ( int d=0; d<12; d++ ){
      arg4_l[d] = ZERO_double;
    }
    double arg5_l[12];
    for ( int d=0; d<12; d++ ){
      arg5_l[d] = ZERO_double;
    }
    int map2idx;
    int map3idx;
    map2idx = opDat2Map[n + set_size * 0];
    map3idx = opDat2Map[n + set_size * 1];
    const double* arg2_vec[] = {
       &ind_arg0[10 * map2idx],
       &ind_arg0[10 * map3idx]};
    double* arg4_vec[] = {
      arg4_l,
      arg5_l};

    //user-supplied kernel call
    ls_advec_edges_gpu(arg0+n*2,
                   arg1+n*1,
                   arg2_vec,
                   arg4_vec);
    atomicAdd(&ind_arg1[0+map2idx*12],arg4_l[0]);
    atomicAdd(&ind_arg1[1+map2idx*12],arg4_l[1]);
    atomicAdd(&ind_arg1[2+map2idx*12],arg4_l[2]);
    atomicAdd(&ind_arg1[3+map2idx*12],arg4_l[3]);
    atomicAdd(&ind_arg1[4+map2idx*12],arg4_l[4]);
    atomicAdd(&ind_arg1[5+map2idx*12],arg4_l[5]);
    atomicAdd(&ind_arg1[6+map2idx*12],arg4_l[6]);
    atomicAdd(&ind_arg1[7+map2idx*12],arg4_l[7]);
    atomicAdd(&ind_arg1[8+map2idx*12],arg4_l[8]);
    atomicAdd(&ind_arg1[9+map2idx*12],arg4_l[9]);
    atomicAdd(&ind_arg1[10+map2idx*12],arg4_l[10]);
    atomicAdd(&ind_arg1[11+map2idx*12],arg4_l[11]);
    atomicAdd(&ind_arg1[0+map3idx*12],arg5_l[0]);
    atomicAdd(&ind_arg1[1+map3idx*12],arg5_l[1]);
    atomicAdd(&ind_arg1[2+map3idx*12],arg5_l[2]);
    atomicAdd(&ind_arg1[3+map3idx*12],arg5_l[3]);
    atomicAdd(&ind_arg1[4+map3idx*12],arg5_l[4]);
    atomicAdd(&ind_arg1[5+map3idx*12],arg5_l[5]);
    atomicAdd(&ind_arg1[6+map3idx*12],arg5_l[6]);
    atomicAdd(&ind_arg1[7+map3idx*12],arg5_l[7]);
    atomicAdd(&ind_arg1[8+map3idx*12],arg5_l[8]);
    atomicAdd(&ind_arg1[9+map3idx*12],arg5_l[9]);
    atomicAdd(&ind_arg1[10+map3idx*12],arg5_l[10]);
    atomicAdd(&ind_arg1[11+map3idx*12],arg5_l[11]);
  }
}


//host stub function
void op_par_loop_ls_advec_edges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 10, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 12, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(53);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[53].name      = name;
  OP_kernels[53].count    += 1;


  int    ninds   = 2;
  int    inds[6] = {-1,-1,0,0,1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_advec_edges\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_53
      int nthread = OP_BLOCK_SIZE_53;
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
        op_cuda_ls_advec_edges<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg4.data_d,
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
  OP_kernels[53].time     += wall_t2 - wall_t1;
}