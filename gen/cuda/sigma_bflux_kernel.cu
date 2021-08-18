//
// auto-generated by op2.py
//

//user function
__device__ void sigma_bflux_gpu( const int *bedgeNum, const double *sJ, const double *nx,
                        const double *ny, const double *s, double *sigFx,
                        double *sigFy) {
  int exInd = *bedgeNum * 6;

  for(int i = 0; i < 6; i++) {
    sigFx[exInd + i] += gaussW_g_cuda[i] * sJ[exInd + i] * nx[exInd + i] * s[exInd + i];
    sigFy[exInd + i] += gaussW_g_cuda[i] * sJ[exInd + i] * ny[exInd + i] * s[exInd + i];
  }

}

// CUDA kernel function
__global__ void op_cuda_sigma_bflux(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  double *__restrict ind_arg4,
  double *__restrict ind_arg5,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg5_l[18];
    for ( int d=0; d<18; d++ ){
      arg5_l[d] = ZERO_double;
    }
    double arg6_l[18];
    for ( int d=0; d<18; d++ ){
      arg6_l[d] = ZERO_double;
    }
    int map1idx;
    map1idx = opDat1Map[n + set_size * 0];

    //user-supplied kernel call
    sigma_bflux_gpu(arg0+n*1,
                ind_arg0+map1idx*18,
                ind_arg1+map1idx*18,
                ind_arg2+map1idx*18,
                ind_arg3+map1idx*18,
                arg5_l,
                arg6_l);
    atomicAdd(&ind_arg4[0+map1idx*18],arg5_l[0]);
    atomicAdd(&ind_arg4[1+map1idx*18],arg5_l[1]);
    atomicAdd(&ind_arg4[2+map1idx*18],arg5_l[2]);
    atomicAdd(&ind_arg4[3+map1idx*18],arg5_l[3]);
    atomicAdd(&ind_arg4[4+map1idx*18],arg5_l[4]);
    atomicAdd(&ind_arg4[5+map1idx*18],arg5_l[5]);
    atomicAdd(&ind_arg4[6+map1idx*18],arg5_l[6]);
    atomicAdd(&ind_arg4[7+map1idx*18],arg5_l[7]);
    atomicAdd(&ind_arg4[8+map1idx*18],arg5_l[8]);
    atomicAdd(&ind_arg4[9+map1idx*18],arg5_l[9]);
    atomicAdd(&ind_arg4[10+map1idx*18],arg5_l[10]);
    atomicAdd(&ind_arg4[11+map1idx*18],arg5_l[11]);
    atomicAdd(&ind_arg4[12+map1idx*18],arg5_l[12]);
    atomicAdd(&ind_arg4[13+map1idx*18],arg5_l[13]);
    atomicAdd(&ind_arg4[14+map1idx*18],arg5_l[14]);
    atomicAdd(&ind_arg4[15+map1idx*18],arg5_l[15]);
    atomicAdd(&ind_arg4[16+map1idx*18],arg5_l[16]);
    atomicAdd(&ind_arg4[17+map1idx*18],arg5_l[17]);
    atomicAdd(&ind_arg5[0+map1idx*18],arg6_l[0]);
    atomicAdd(&ind_arg5[1+map1idx*18],arg6_l[1]);
    atomicAdd(&ind_arg5[2+map1idx*18],arg6_l[2]);
    atomicAdd(&ind_arg5[3+map1idx*18],arg6_l[3]);
    atomicAdd(&ind_arg5[4+map1idx*18],arg6_l[4]);
    atomicAdd(&ind_arg5[5+map1idx*18],arg6_l[5]);
    atomicAdd(&ind_arg5[6+map1idx*18],arg6_l[6]);
    atomicAdd(&ind_arg5[7+map1idx*18],arg6_l[7]);
    atomicAdd(&ind_arg5[8+map1idx*18],arg6_l[8]);
    atomicAdd(&ind_arg5[9+map1idx*18],arg6_l[9]);
    atomicAdd(&ind_arg5[10+map1idx*18],arg6_l[10]);
    atomicAdd(&ind_arg5[11+map1idx*18],arg6_l[11]);
    atomicAdd(&ind_arg5[12+map1idx*18],arg6_l[12]);
    atomicAdd(&ind_arg5[13+map1idx*18],arg6_l[13]);
    atomicAdd(&ind_arg5[14+map1idx*18],arg6_l[14]);
    atomicAdd(&ind_arg5[15+map1idx*18],arg6_l[15]);
    atomicAdd(&ind_arg5[16+map1idx*18],arg6_l[16]);
    atomicAdd(&ind_arg5[17+map1idx*18],arg6_l[17]);
  }
}


//host stub function
void op_par_loop_sigma_bflux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6){

  int nargs = 7;
  op_arg args[7];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(60);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[60].name      = name;
  OP_kernels[60].count    += 1;


  int    ninds   = 6;
  int    inds[7] = {-1,0,1,2,3,4,5};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: sigma_bflux\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_60
      int nthread = OP_BLOCK_SIZE_60;
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
        op_cuda_sigma_bflux<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg2.data_d,
        (double *)arg3.data_d,
        (double *)arg4.data_d,
        (double *)arg5.data_d,
        (double *)arg6.data_d,
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
  OP_kernels[60].time     += wall_t2 - wall_t1;
}
