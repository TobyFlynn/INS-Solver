//
// auto-generated by op2.py
//

//user function
__device__ void poisson_edges_gpu( const double *uL, const double *opL, double *rhsL,
                          const double *uR, const double *opR, double *rhsR) {
  for(int m = 0; m < 15; m++) {
    int ind = m * 15;
    for(int n = 0; n < 15; n++) {
      rhsL[m] += opL[ind + n] * uR[n];
      rhsR[m] += opR[ind + n] * uL[n];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_edges(
  const double *__restrict ind_arg0,
  double *__restrict ind_arg1,
  const int *__restrict opDat0Map,
  const double *__restrict arg1,
  const double *__restrict arg4,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg2_l[15];
    for ( int d=0; d<15; d++ ){
      arg2_l[d] = ZERO_double;
    }
    double arg5_l[15];
    for ( int d=0; d<15; d++ ){
      arg5_l[d] = ZERO_double;
    }
    int map0idx;
    int map3idx;
    map0idx = opDat0Map[n + set_size * 0];
    map3idx = opDat0Map[n + set_size * 1];

    //user-supplied kernel call
    poisson_edges_gpu(ind_arg0+map0idx*15,
                  arg1+n*225,
                  arg2_l,
                  ind_arg0+map3idx*15,
                  arg4+n*225,
                  arg5_l);
    atomicAdd(&ind_arg1[0+map0idx*15],arg2_l[0]);
    atomicAdd(&ind_arg1[1+map0idx*15],arg2_l[1]);
    atomicAdd(&ind_arg1[2+map0idx*15],arg2_l[2]);
    atomicAdd(&ind_arg1[3+map0idx*15],arg2_l[3]);
    atomicAdd(&ind_arg1[4+map0idx*15],arg2_l[4]);
    atomicAdd(&ind_arg1[5+map0idx*15],arg2_l[5]);
    atomicAdd(&ind_arg1[6+map0idx*15],arg2_l[6]);
    atomicAdd(&ind_arg1[7+map0idx*15],arg2_l[7]);
    atomicAdd(&ind_arg1[8+map0idx*15],arg2_l[8]);
    atomicAdd(&ind_arg1[9+map0idx*15],arg2_l[9]);
    atomicAdd(&ind_arg1[10+map0idx*15],arg2_l[10]);
    atomicAdd(&ind_arg1[11+map0idx*15],arg2_l[11]);
    atomicAdd(&ind_arg1[12+map0idx*15],arg2_l[12]);
    atomicAdd(&ind_arg1[13+map0idx*15],arg2_l[13]);
    atomicAdd(&ind_arg1[14+map0idx*15],arg2_l[14]);
    atomicAdd(&ind_arg1[0+map3idx*15],arg5_l[0]);
    atomicAdd(&ind_arg1[1+map3idx*15],arg5_l[1]);
    atomicAdd(&ind_arg1[2+map3idx*15],arg5_l[2]);
    atomicAdd(&ind_arg1[3+map3idx*15],arg5_l[3]);
    atomicAdd(&ind_arg1[4+map3idx*15],arg5_l[4]);
    atomicAdd(&ind_arg1[5+map3idx*15],arg5_l[5]);
    atomicAdd(&ind_arg1[6+map3idx*15],arg5_l[6]);
    atomicAdd(&ind_arg1[7+map3idx*15],arg5_l[7]);
    atomicAdd(&ind_arg1[8+map3idx*15],arg5_l[8]);
    atomicAdd(&ind_arg1[9+map3idx*15],arg5_l[9]);
    atomicAdd(&ind_arg1[10+map3idx*15],arg5_l[10]);
    atomicAdd(&ind_arg1[11+map3idx*15],arg5_l[11]);
    atomicAdd(&ind_arg1[12+map3idx*15],arg5_l[12]);
    atomicAdd(&ind_arg1[13+map3idx*15],arg5_l[13]);
    atomicAdd(&ind_arg1[14+map3idx*15],arg5_l[14]);
  }
}


//host stub function
void op_par_loop_poisson_edges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(17);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[17].name      = name;
  OP_kernels[17].count    += 1;


  int    ninds   = 2;
  int    inds[6] = {0,-1,1,0,-1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_edges\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_17
      int nthread = OP_BLOCK_SIZE_17;
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
        op_cuda_poisson_edges<<<nblocks,nthread>>>(
        (double *)arg0.data_d,
        (double *)arg2.data_d,
        arg0.map_data_d,
        (double*)arg1.data_d,
        (double*)arg4.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[17].time     += wall_t2 - wall_t1;
}
