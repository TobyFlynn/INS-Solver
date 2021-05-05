//
// auto-generated by op2.py
//

//user function
__device__ void tau_gpu( const int *edgeNum, const double **J, const double **sJ, double **tau) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  double hminL = 2.0 * J[0][FMASK_cuda[edgeL * 5]] / sJ[0][edgeL * 5];
  double hminR = 2.0 * J[1][FMASK_cuda[edgeR * 5]] / sJ[1][edgeR * 5];

  if(hminL > hminR) {
    tau[0][edgeL] += 100.0 * 15.0 / hminR;
    tau[1][edgeR] += 100.0 * 15.0 / hminR;
  } else {
    tau[0][edgeL] += 100.0 * 15.0 / hminL;
    tau[1][edgeR] += 100.0 * 15.0 / hminL;
  }

}

// CUDA kernel function
__global__ void op_cuda_tau(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  double *__restrict ind_arg2,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  int start,
  int end,
  int   set_size) {
  double arg5_l[3];
  double arg6_l[3];
  double *arg5_vec[2] = {
    arg5_l,
    arg6_l,
  };
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg5_l[3];
    for ( int d=0; d<3; d++ ){
      arg5_l[d] = ZERO_double;
    }
    double arg6_l[3];
    for ( int d=0; d<3; d++ ){
      arg6_l[d] = ZERO_double;
    }
    int map1idx;
    int map2idx;
    map1idx = opDat1Map[n + set_size * 0];
    map2idx = opDat1Map[n + set_size * 1];
    const double* arg1_vec[] = {
       &ind_arg0[15 * map1idx],
       &ind_arg0[15 * map2idx]};
    const double* arg3_vec[] = {
       &ind_arg1[15 * map1idx],
       &ind_arg1[15 * map2idx]};
    double* arg5_vec[] = {
       &ind_arg2[3 * map1idx],
       &ind_arg2[3 * map2idx]};

    //user-supplied kernel call
    tau_gpu(arg0+n*2,
        arg1_vec,
        arg3_vec,
        arg5_vec);
    atomicAdd(&ind_arg2[0+map1idx*3],arg5_l[0]);
    atomicAdd(&ind_arg2[1+map1idx*3],arg5_l[1]);
    atomicAdd(&ind_arg2[2+map1idx*3],arg5_l[2]);
    atomicAdd(&ind_arg2[0+map2idx*3],arg6_l[0]);
    atomicAdd(&ind_arg2[1+map2idx*3],arg6_l[1]);
    atomicAdd(&ind_arg2[2+map2idx*3],arg6_l[2]);
  }
}


//host stub function
void op_par_loop_tau(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3,
  op_arg arg5){

  int nargs = 7;
  op_arg args[7];

  args[0] = arg0;
  arg1.idx = 0;
  args[1] = arg1;
  for ( int v=1; v<2; v++ ){
    args[1 + v] = op_arg_dat(arg1.dat, v, arg1.map, 15, "double", OP_READ);
  }

  arg3.idx = 0;
  args[3] = arg3;
  for ( int v=1; v<2; v++ ){
    args[3 + v] = op_arg_dat(arg3.dat, v, arg3.map, 15, "double", OP_READ);
  }

  arg5.idx = 0;
  args[5] = arg5;
  for ( int v=1; v<2; v++ ){
    args[5 + v] = op_arg_dat(arg5.dat, v, arg5.map, 3, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(21);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[21].name      = name;
  OP_kernels[21].count    += 1;


  int    ninds   = 3;
  int    inds[7] = {-1,0,0,1,1,2,2};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: tau\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_21
      int nthread = OP_BLOCK_SIZE_21;
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
        op_cuda_tau<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg3.data_d,
        (double *)arg5.data_d,
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
  OP_kernels[21].time     += wall_t2 - wall_t1;
}
