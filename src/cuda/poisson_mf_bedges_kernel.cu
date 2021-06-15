//
// auto-generated by op2.py
//

//user function
__device__ void poisson_mf_bedges_gpu( const int *bedgeType, const int *bedgeNum,
                              const int *d0, const int *d1, const int *d2,
                              const double *sJ, const double *nx,
                              const double *ny, const double *tau,
                              const double *u, const double *dudx,
                              const double *dudy, double *fluxX, double *fluxY,
                              double *flux) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 2 * 7;

  if(*bedgeType == *d0 || *bedgeType == *d1 || *bedgeType == *d2) {
    for(int i = 0; i < 7; i++) {
      flux[exInd + i] += gaussW_g_cuda[i] * sJ[exInd + i] * (nx[exInd + i] * dudx[exInd + i] + ny[exInd + i] * dudy[exInd + i] - tau[*bedgeNum] * (u[exInd + i]));
    }
  } else {
    for(int i = 0; i < 7; i++) {
      fluxX[exInd + i] += nx[exInd + i] * gaussW_g_cuda[i] * sJ[exInd + i] * u[exInd + i];
      fluxY[exInd + i] += ny[exInd + i] * gaussW_g_cuda[i] * sJ[exInd + i] * u[exInd + i];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_mf_bedges(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  const double *__restrict ind_arg6,
  double *__restrict ind_arg7,
  double *__restrict ind_arg8,
  double *__restrict ind_arg9,
  const int *__restrict opDat5Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  const int *arg2,
  const int *arg3,
  const int *arg4,
  int start,
  int end,
  int   set_size) {
  double arg12_l[21];
  double arg13_l[21];
  double arg14_l[21];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg12_l[21];
    for ( int d=0; d<21; d++ ){
      arg12_l[d] = ZERO_double;
    }
    double arg13_l[21];
    for ( int d=0; d<21; d++ ){
      arg13_l[d] = ZERO_double;
    }
    double arg14_l[21];
    for ( int d=0; d<21; d++ ){
      arg14_l[d] = ZERO_double;
    }
    int map5idx;
    map5idx = opDat5Map[n + set_size * 0];

    //user-supplied kernel call
    poisson_mf_bedges_gpu(arg0+n*1,
                      arg1+n*1,
                      arg2,
                      arg3,
                      arg4,
                      ind_arg0+map5idx*21,
                      ind_arg1+map5idx*21,
                      ind_arg2+map5idx*21,
                      ind_arg3+map5idx*3,
                      ind_arg4+map5idx*21,
                      ind_arg5+map5idx*21,
                      ind_arg6+map5idx*21,
                      arg12_l,
                      arg13_l,
                      arg14_l);
    atomicAdd(&ind_arg7[0+map5idx*21],arg12_l[0]);
    atomicAdd(&ind_arg7[1+map5idx*21],arg12_l[1]);
    atomicAdd(&ind_arg7[2+map5idx*21],arg12_l[2]);
    atomicAdd(&ind_arg7[3+map5idx*21],arg12_l[3]);
    atomicAdd(&ind_arg7[4+map5idx*21],arg12_l[4]);
    atomicAdd(&ind_arg7[5+map5idx*21],arg12_l[5]);
    atomicAdd(&ind_arg7[6+map5idx*21],arg12_l[6]);
    atomicAdd(&ind_arg7[7+map5idx*21],arg12_l[7]);
    atomicAdd(&ind_arg7[8+map5idx*21],arg12_l[8]);
    atomicAdd(&ind_arg7[9+map5idx*21],arg12_l[9]);
    atomicAdd(&ind_arg7[10+map5idx*21],arg12_l[10]);
    atomicAdd(&ind_arg7[11+map5idx*21],arg12_l[11]);
    atomicAdd(&ind_arg7[12+map5idx*21],arg12_l[12]);
    atomicAdd(&ind_arg7[13+map5idx*21],arg12_l[13]);
    atomicAdd(&ind_arg7[14+map5idx*21],arg12_l[14]);
    atomicAdd(&ind_arg7[15+map5idx*21],arg12_l[15]);
    atomicAdd(&ind_arg7[16+map5idx*21],arg12_l[16]);
    atomicAdd(&ind_arg7[17+map5idx*21],arg12_l[17]);
    atomicAdd(&ind_arg7[18+map5idx*21],arg12_l[18]);
    atomicAdd(&ind_arg7[19+map5idx*21],arg12_l[19]);
    atomicAdd(&ind_arg7[20+map5idx*21],arg12_l[20]);
    atomicAdd(&ind_arg8[0+map5idx*21],arg13_l[0]);
    atomicAdd(&ind_arg8[1+map5idx*21],arg13_l[1]);
    atomicAdd(&ind_arg8[2+map5idx*21],arg13_l[2]);
    atomicAdd(&ind_arg8[3+map5idx*21],arg13_l[3]);
    atomicAdd(&ind_arg8[4+map5idx*21],arg13_l[4]);
    atomicAdd(&ind_arg8[5+map5idx*21],arg13_l[5]);
    atomicAdd(&ind_arg8[6+map5idx*21],arg13_l[6]);
    atomicAdd(&ind_arg8[7+map5idx*21],arg13_l[7]);
    atomicAdd(&ind_arg8[8+map5idx*21],arg13_l[8]);
    atomicAdd(&ind_arg8[9+map5idx*21],arg13_l[9]);
    atomicAdd(&ind_arg8[10+map5idx*21],arg13_l[10]);
    atomicAdd(&ind_arg8[11+map5idx*21],arg13_l[11]);
    atomicAdd(&ind_arg8[12+map5idx*21],arg13_l[12]);
    atomicAdd(&ind_arg8[13+map5idx*21],arg13_l[13]);
    atomicAdd(&ind_arg8[14+map5idx*21],arg13_l[14]);
    atomicAdd(&ind_arg8[15+map5idx*21],arg13_l[15]);
    atomicAdd(&ind_arg8[16+map5idx*21],arg13_l[16]);
    atomicAdd(&ind_arg8[17+map5idx*21],arg13_l[17]);
    atomicAdd(&ind_arg8[18+map5idx*21],arg13_l[18]);
    atomicAdd(&ind_arg8[19+map5idx*21],arg13_l[19]);
    atomicAdd(&ind_arg8[20+map5idx*21],arg13_l[20]);
    atomicAdd(&ind_arg9[0+map5idx*21],arg14_l[0]);
    atomicAdd(&ind_arg9[1+map5idx*21],arg14_l[1]);
    atomicAdd(&ind_arg9[2+map5idx*21],arg14_l[2]);
    atomicAdd(&ind_arg9[3+map5idx*21],arg14_l[3]);
    atomicAdd(&ind_arg9[4+map5idx*21],arg14_l[4]);
    atomicAdd(&ind_arg9[5+map5idx*21],arg14_l[5]);
    atomicAdd(&ind_arg9[6+map5idx*21],arg14_l[6]);
    atomicAdd(&ind_arg9[7+map5idx*21],arg14_l[7]);
    atomicAdd(&ind_arg9[8+map5idx*21],arg14_l[8]);
    atomicAdd(&ind_arg9[9+map5idx*21],arg14_l[9]);
    atomicAdd(&ind_arg9[10+map5idx*21],arg14_l[10]);
    atomicAdd(&ind_arg9[11+map5idx*21],arg14_l[11]);
    atomicAdd(&ind_arg9[12+map5idx*21],arg14_l[12]);
    atomicAdd(&ind_arg9[13+map5idx*21],arg14_l[13]);
    atomicAdd(&ind_arg9[14+map5idx*21],arg14_l[14]);
    atomicAdd(&ind_arg9[15+map5idx*21],arg14_l[15]);
    atomicAdd(&ind_arg9[16+map5idx*21],arg14_l[16]);
    atomicAdd(&ind_arg9[17+map5idx*21],arg14_l[17]);
    atomicAdd(&ind_arg9[18+map5idx*21],arg14_l[18]);
    atomicAdd(&ind_arg9[19+map5idx*21],arg14_l[19]);
    atomicAdd(&ind_arg9[20+map5idx*21],arg14_l[20]);
  }
}


//host stub function
void op_par_loop_poisson_mf_bedges(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14){

  int*arg2h = (int *)arg2.data;
  int*arg3h = (int *)arg3.data;
  int*arg4h = (int *)arg4.data;
  int nargs = 15;
  op_arg args[15];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  args[10] = arg10;
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(34);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[34].name      = name;
  OP_kernels[34].count    += 1;


  int    ninds   = 10;
  int    inds[15] = {-1,-1,-1,-1,-1,0,1,2,3,4,5,6,7,8,9};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_mf_bedges\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(int));
    consts_bytes += ROUND_UP(1*sizeof(int));
    consts_bytes += ROUND_UP(1*sizeof(int));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg2.data   = OP_consts_h + consts_bytes;
    arg2.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg2.data)[d] = arg2h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    arg3.data   = OP_consts_h + consts_bytes;
    arg3.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg3.data)[d] = arg3h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    arg4.data   = OP_consts_h + consts_bytes;
    arg4.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg4.data)[d] = arg4h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_34
      int nthread = OP_BLOCK_SIZE_34;
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
        op_cuda_poisson_mf_bedges<<<nblocks,nthread>>>(
        (double *)arg5.data_d,
        (double *)arg6.data_d,
        (double *)arg7.data_d,
        (double *)arg8.data_d,
        (double *)arg9.data_d,
        (double *)arg10.data_d,
        (double *)arg11.data_d,
        (double *)arg12.data_d,
        (double *)arg13.data_d,
        (double *)arg14.data_d,
        arg5.map_data_d,
        (int*)arg0.data_d,
        (int*)arg1.data_d,
        (int*)arg2.data_d,
        (int*)arg3.data_d,
        (int*)arg4.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[34].time     += wall_t2 - wall_t1;
}
