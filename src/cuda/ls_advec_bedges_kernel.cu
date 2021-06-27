//
// auto-generated by op2.py
//

//user function
__device__ void ls_advec_bedges_gpu( const int *bedge_type, const int *bedgeNum,
                            const double *x, const double *y, const double *q,
                            double *exQ) {
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
    exQ[exInd + i] += q[fmask[i]];
  }

}

// CUDA kernel function
__global__ void op_cuda_ls_advec_bedges(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  double *__restrict ind_arg3,
  const int *__restrict opDat2Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  int start,
  int end,
  int   set_size) {
  double arg5_l[15];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg5_l[15];
    for ( int d=0; d<15; d++ ){
      arg5_l[d] = ZERO_double;
    }
    int map2idx;
    map2idx = opDat2Map[n + set_size * 0];

    //user-supplied kernel call
    ls_advec_bedges_gpu(arg0+n*1,
                    arg1+n*1,
                    ind_arg0+map2idx*15,
                    ind_arg1+map2idx*15,
                    ind_arg2+map2idx*15,
                    arg5_l);
    atomicAdd(&ind_arg3[0+map2idx*15],arg5_l[0]);
    atomicAdd(&ind_arg3[1+map2idx*15],arg5_l[1]);
    atomicAdd(&ind_arg3[2+map2idx*15],arg5_l[2]);
    atomicAdd(&ind_arg3[3+map2idx*15],arg5_l[3]);
    atomicAdd(&ind_arg3[4+map2idx*15],arg5_l[4]);
    atomicAdd(&ind_arg3[5+map2idx*15],arg5_l[5]);
    atomicAdd(&ind_arg3[6+map2idx*15],arg5_l[6]);
    atomicAdd(&ind_arg3[7+map2idx*15],arg5_l[7]);
    atomicAdd(&ind_arg3[8+map2idx*15],arg5_l[8]);
    atomicAdd(&ind_arg3[9+map2idx*15],arg5_l[9]);
    atomicAdd(&ind_arg3[10+map2idx*15],arg5_l[10]);
    atomicAdd(&ind_arg3[11+map2idx*15],arg5_l[11]);
    atomicAdd(&ind_arg3[12+map2idx*15],arg5_l[12]);
    atomicAdd(&ind_arg3[13+map2idx*15],arg5_l[13]);
    atomicAdd(&ind_arg3[14+map2idx*15],arg5_l[14]);
  }
}


//host stub function
void op_par_loop_ls_advec_bedges(char const *name, op_set set,
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
  op_timing_realloc(57);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[57].name      = name;
  OP_kernels[57].count    += 1;


  int    ninds   = 4;
  int    inds[6] = {-1,-1,0,1,2,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: ls_advec_bedges\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_57
      int nthread = OP_BLOCK_SIZE_57;
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
        op_cuda_ls_advec_bedges<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg3.data_d,
        (double *)arg4.data_d,
        (double *)arg5.data_d,
        arg2.map_data_d,
        (int*)arg0.data_d,
        (int*)arg1.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[57].time     += wall_t2 - wall_t1;
}
