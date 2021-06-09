//
// auto-generated by op2.py
//

//user function
__device__ void init_edges_gpu( const int *edgeNum, const double **x, const double **y,
                       bool *reverse) {
  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  if(edgeR == 0) {
    if(edgeL == 0) {
      *reverse = !(x[0][0] == x[1][0] && y[0][0] == y[1][0]);
    } else if(edgeL == 1) {
      *reverse = !(x[0][1] == x[1][0] && y[0][1] == y[1][0]);
    } else {
      *reverse = !(x[0][2] == x[1][0] && y[0][2] == y[1][0]);
    }
  } else if(edgeR == 1) {
    if(edgeL == 0) {
      *reverse = !(x[0][0] == x[1][1] && y[0][0] == y[1][1]);
    } else if(edgeL == 1) {
      *reverse = !(x[0][1] == x[1][1] && y[0][1] == y[1][1]);
    } else {
      *reverse = !(x[0][2] == x[1][1] && y[0][2] == y[1][1]);
    }
  } else {
    if(edgeL == 0) {
      *reverse = !(x[0][0] == x[1][2] && y[0][0] == y[1][2]);
    } else if(edgeL == 1) {
      *reverse = !(x[0][1] == x[1][2] && y[0][1] == y[1][2]);
    } else {
      *reverse = !(x[0][2] == x[1][2] && y[0][2] == y[1][2]);
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_init_edges(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  bool *arg5,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
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

    //user-supplied kernel call
    init_edges_gpu(arg0+n*2,
               arg1_vec,
               arg3_vec,
               arg5+n*1);
  }
}


//host stub function
void op_par_loop_init_edges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg3,
  op_arg arg5){

  int nargs = 6;
  op_arg args[6];

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

  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(2);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[2].name      = name;
  OP_kernels[2].count    += 1;


  int    ninds   = 2;
  int    inds[6] = {-1,0,0,1,1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: init_edges\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_2
      int nthread = OP_BLOCK_SIZE_2;
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
        op_cuda_init_edges<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg3.data_d,
        arg1.map_data_d,
        (int*)arg0.data_d,
        (bool*)arg5.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[2].time     += wall_t2 - wall_t1;
}