//
// auto-generated by op2.py
//

//user function
__device__ void advection_faces_gpu( const int *edgeNum, const bool *rev, const double **q0,
                            const double **q1, double **exQ0, double **exQ1) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

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
    exQ0[0][exInd + i] += q0[1][rInd];
    exQ1[0][exInd + i] += q1[1][rInd];
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
    exQ0[1][exInd + i] += q0[0][lInd];
    exQ1[1][exInd + i] += q1[0][lInd];
  }

}

// CUDA kernel function
__global__ void op_cuda_advection_faces(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  double *__restrict ind_arg2,
  double *__restrict ind_arg3,
  const int *__restrict opDat2Map,
  const int *__restrict arg0,
  const bool *__restrict arg1,
  int start,
  int end,
  int   set_size) {
  double arg6_l[15];
  double arg7_l[15];
  double arg8_l[15];
  double arg9_l[15];
  double *arg6_vec[2] = {
    arg6_l,
    arg7_l,
  };
  double *arg8_vec[2] = {
    arg8_l,
    arg9_l,
  };
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg6_l[15];
    for ( int d=0; d<15; d++ ){
      arg6_l[d] = ZERO_double;
    }
    double arg7_l[15];
    for ( int d=0; d<15; d++ ){
      arg7_l[d] = ZERO_double;
    }
    double arg8_l[15];
    for ( int d=0; d<15; d++ ){
      arg8_l[d] = ZERO_double;
    }
    double arg9_l[15];
    for ( int d=0; d<15; d++ ){
      arg9_l[d] = ZERO_double;
    }
    int map2idx;
    int map3idx;
    map2idx = opDat2Map[n + set_size * 0];
    map3idx = opDat2Map[n + set_size * 1];
    const double* arg2_vec[] = {
       &ind_arg0[15 * map2idx],
       &ind_arg0[15 * map3idx]};
    const double* arg4_vec[] = {
       &ind_arg1[15 * map2idx],
       &ind_arg1[15 * map3idx]};
    double* arg6_vec[] = {
       &ind_arg2[15 * map2idx],
       &ind_arg2[15 * map3idx]};
    double* arg8_vec[] = {
       &ind_arg3[15 * map2idx],
       &ind_arg3[15 * map3idx]};

    //user-supplied kernel call
    advection_faces_gpu(arg0+n*2,
                    arg1+n*1,
                    arg2_vec,
                    arg4_vec,
                    arg6_vec,
                    arg8_vec);
    atomicAdd(&ind_arg2[0+map2idx*15],arg6_l[0]);
    atomicAdd(&ind_arg2[1+map2idx*15],arg6_l[1]);
    atomicAdd(&ind_arg2[2+map2idx*15],arg6_l[2]);
    atomicAdd(&ind_arg2[3+map2idx*15],arg6_l[3]);
    atomicAdd(&ind_arg2[4+map2idx*15],arg6_l[4]);
    atomicAdd(&ind_arg2[5+map2idx*15],arg6_l[5]);
    atomicAdd(&ind_arg2[6+map2idx*15],arg6_l[6]);
    atomicAdd(&ind_arg2[7+map2idx*15],arg6_l[7]);
    atomicAdd(&ind_arg2[8+map2idx*15],arg6_l[8]);
    atomicAdd(&ind_arg2[9+map2idx*15],arg6_l[9]);
    atomicAdd(&ind_arg2[10+map2idx*15],arg6_l[10]);
    atomicAdd(&ind_arg2[11+map2idx*15],arg6_l[11]);
    atomicAdd(&ind_arg2[12+map2idx*15],arg6_l[12]);
    atomicAdd(&ind_arg2[13+map2idx*15],arg6_l[13]);
    atomicAdd(&ind_arg2[14+map2idx*15],arg6_l[14]);
    atomicAdd(&ind_arg2[0+map3idx*15],arg7_l[0]);
    atomicAdd(&ind_arg2[1+map3idx*15],arg7_l[1]);
    atomicAdd(&ind_arg2[2+map3idx*15],arg7_l[2]);
    atomicAdd(&ind_arg2[3+map3idx*15],arg7_l[3]);
    atomicAdd(&ind_arg2[4+map3idx*15],arg7_l[4]);
    atomicAdd(&ind_arg2[5+map3idx*15],arg7_l[5]);
    atomicAdd(&ind_arg2[6+map3idx*15],arg7_l[6]);
    atomicAdd(&ind_arg2[7+map3idx*15],arg7_l[7]);
    atomicAdd(&ind_arg2[8+map3idx*15],arg7_l[8]);
    atomicAdd(&ind_arg2[9+map3idx*15],arg7_l[9]);
    atomicAdd(&ind_arg2[10+map3idx*15],arg7_l[10]);
    atomicAdd(&ind_arg2[11+map3idx*15],arg7_l[11]);
    atomicAdd(&ind_arg2[12+map3idx*15],arg7_l[12]);
    atomicAdd(&ind_arg2[13+map3idx*15],arg7_l[13]);
    atomicAdd(&ind_arg2[14+map3idx*15],arg7_l[14]);
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
    atomicAdd(&ind_arg3[0+map3idx*15],arg9_l[0]);
    atomicAdd(&ind_arg3[1+map3idx*15],arg9_l[1]);
    atomicAdd(&ind_arg3[2+map3idx*15],arg9_l[2]);
    atomicAdd(&ind_arg3[3+map3idx*15],arg9_l[3]);
    atomicAdd(&ind_arg3[4+map3idx*15],arg9_l[4]);
    atomicAdd(&ind_arg3[5+map3idx*15],arg9_l[5]);
    atomicAdd(&ind_arg3[6+map3idx*15],arg9_l[6]);
    atomicAdd(&ind_arg3[7+map3idx*15],arg9_l[7]);
    atomicAdd(&ind_arg3[8+map3idx*15],arg9_l[8]);
    atomicAdd(&ind_arg3[9+map3idx*15],arg9_l[9]);
    atomicAdd(&ind_arg3[10+map3idx*15],arg9_l[10]);
    atomicAdd(&ind_arg3[11+map3idx*15],arg9_l[11]);
    atomicAdd(&ind_arg3[12+map3idx*15],arg9_l[12]);
    atomicAdd(&ind_arg3[13+map3idx*15],arg9_l[13]);
    atomicAdd(&ind_arg3[14+map3idx*15],arg9_l[14]);
  }
}


//host stub function
void op_par_loop_advection_faces(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg4,
  op_arg arg6,
  op_arg arg8){

  int nargs = 10;
  op_arg args[10];

  args[0] = arg0;
  args[1] = arg1;
  arg2.idx = 0;
  args[2] = arg2;
  for ( int v=1; v<2; v++ ){
    args[2 + v] = op_arg_dat(arg2.dat, v, arg2.map, 15, "double", OP_READ);
  }

  arg4.idx = 0;
  args[4] = arg4;
  for ( int v=1; v<2; v++ ){
    args[4 + v] = op_arg_dat(arg4.dat, v, arg4.map, 15, "double", OP_READ);
  }

  arg6.idx = 0;
  args[6] = arg6;
  for ( int v=1; v<2; v++ ){
    args[6 + v] = op_arg_dat(arg6.dat, v, arg6.map, 15, "double", OP_INC);
  }

  arg8.idx = 0;
  args[8] = arg8;
  for ( int v=1; v<2; v++ ){
    args[8 + v] = op_arg_dat(arg8.dat, v, arg8.map, 15, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(49);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[49].name      = name;
  OP_kernels[49].count    += 1;


  int    ninds   = 4;
  int    inds[10] = {-1,-1,0,0,1,1,2,2,3,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: advection_faces\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_49
      int nthread = OP_BLOCK_SIZE_49;
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
        op_cuda_advection_faces<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg4.data_d,
        (double *)arg6.data_d,
        (double *)arg8.data_d,
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
  OP_kernels[49].time     += wall_t2 - wall_t1;
}
