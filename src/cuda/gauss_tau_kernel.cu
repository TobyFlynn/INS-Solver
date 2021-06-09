//
// auto-generated by op2.py
//

//user function
__device__ void gauss_tau_gpu( const int *edgeNum, const double **fscale, double **tau) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  if(fscale[0][edgeL * 5] > fscale[1][edgeR * 5]) {
    tau[0][edgeL] += 20 * 25 * fscale[0][edgeL * 5];
    tau[1][edgeR] += 20 * 25 * fscale[0][edgeL * 5];
  } else {
    tau[0][edgeL] += 20 * 25 * fscale[1][edgeR * 5];
    tau[1][edgeR] += 20 * 25 * fscale[1][edgeR * 5];
  }

}

// CUDA kernel function
__global__ void op_cuda_gauss_tau(
  const double *__restrict ind_arg0,
  double *__restrict ind_arg1,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  int start,
  int end,
  int   set_size) {
  double arg3_l[3];
  double arg4_l[3];
  double *arg3_vec[2] = {
    arg3_l,
    arg4_l,
  };
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg3_l[3];
    for ( int d=0; d<3; d++ ){
      arg3_l[d] = ZERO_double;
    }
    double arg4_l[3];
    for ( int d=0; d<3; d++ ){
      arg4_l[d] = ZERO_double;
    }
    int map1idx;
    int map2idx;
    map1idx = opDat1Map[n + set_size * 0];
    map2idx = opDat1Map[n + set_size * 1];
    const double* arg1_vec[] = {
       &ind_arg0[15 * map1idx],
       &ind_arg0[15 * map2idx]};
    double* arg3_vec[] = {
       &ind_arg1[3 * map1idx],
       &ind_arg1[3 * map2idx]};

    //user-supplied kernel call
    gauss_tau_gpu(arg0+n*2,
              arg1_vec,
              arg3_vec);
    atomicAdd(&ind_arg1[0+map1idx*3],arg3_l[0]);
    atomicAdd(&ind_arg1[1+map1idx*3],arg3_l[1]);
    atomicAdd(&ind_arg1[2+map1idx*3],arg3_l[2]);
    atomicAdd(&ind_arg1[0+map2idx*3],arg4_l[0]);
    atomicAdd(&ind_arg1[1+map2idx*3],arg4_l[1]);
    atomicAdd(&ind_arg1[2+map2idx*3],arg4_l[2]);
  }
}


//host stub function
void op_par_loop_gauss_tau(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  arg1.idx = 0;
  args[1] = arg1;
  for ( int v=1; v<2; v++ ){
    args[1 + v] = op_arg_dat(arg1.dat, v, arg1.map, 15, "double", OP_READ);
  }

  arg3.idx = 0;
  args[3] = arg3;
  for ( int v=1; v<2; v++ ){
    args[3 + v] = op_arg_dat(arg3.dat, v, arg3.map, 3, "double", OP_INC);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(8);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[8].name      = name;
  OP_kernels[8].count    += 1;


  int    ninds   = 2;
  int    inds[5] = {-1,0,0,1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_tau\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_8
      int nthread = OP_BLOCK_SIZE_8;
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
        op_cuda_gauss_tau<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg3.data_d,
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
  OP_kernels[8].time     += wall_t2 - wall_t1;
}
