//
// auto-generated by op2.py
//

//user function
__device__ void set_tau_bc_gpu( const int *bedgeNum, const double *J, const double *sJ,
                       double *tau) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 5;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 5;
  }

  int *fmask;

  if(*bedgeNum == 0) {
    fmask = FMASK_cuda;
  } else if(*bedgeNum == 1) {
    fmask = &FMASK_cuda[5];
  } else {
    fmask = &FMASK_cuda[2 * 5];
  }

  for(int i = 0; i < 5; i++) {
    tau[exInd + i] += 15.0 / (2.0 * J[fmask[i]] / sJ[exInd + i]);
  }

}

// CUDA kernel function
__global__ void op_cuda_set_tau_bc(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  double *__restrict ind_arg2,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  int start,
  int end,
  int   set_size) {
  double arg3_l[15];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg3_l[15];
    for ( int d=0; d<15; d++ ){
      arg3_l[d] = ZERO_double;
    }
    int map1idx;
    map1idx = opDat1Map[n + set_size * 0];

    //user-supplied kernel call
    set_tau_bc_gpu(arg0+n*1,
               ind_arg0+map1idx*15,
               ind_arg1+map1idx*15,
               arg3_l);
    atomicAdd(&ind_arg2[0+map1idx*15],arg3_l[0]);
    atomicAdd(&ind_arg2[1+map1idx*15],arg3_l[1]);
    atomicAdd(&ind_arg2[2+map1idx*15],arg3_l[2]);
    atomicAdd(&ind_arg2[3+map1idx*15],arg3_l[3]);
    atomicAdd(&ind_arg2[4+map1idx*15],arg3_l[4]);
    atomicAdd(&ind_arg2[5+map1idx*15],arg3_l[5]);
    atomicAdd(&ind_arg2[6+map1idx*15],arg3_l[6]);
    atomicAdd(&ind_arg2[7+map1idx*15],arg3_l[7]);
    atomicAdd(&ind_arg2[8+map1idx*15],arg3_l[8]);
    atomicAdd(&ind_arg2[9+map1idx*15],arg3_l[9]);
    atomicAdd(&ind_arg2[10+map1idx*15],arg3_l[10]);
    atomicAdd(&ind_arg2[11+map1idx*15],arg3_l[11]);
    atomicAdd(&ind_arg2[12+map1idx*15],arg3_l[12]);
    atomicAdd(&ind_arg2[13+map1idx*15],arg3_l[13]);
    atomicAdd(&ind_arg2[14+map1idx*15],arg3_l[14]);
  }
}


//host stub function
void op_par_loop_set_tau_bc(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(21);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[21].name      = name;
  OP_kernels[21].count    += 1;


  int    ninds   = 3;
  int    inds[4] = {-1,0,1,2};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: set_tau_bc\n");
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
        op_cuda_set_tau_bc<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg2.data_d,
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
  OP_kernels[21].time     += wall_t2 - wall_t1;
}
