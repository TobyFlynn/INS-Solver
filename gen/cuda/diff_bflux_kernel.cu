//
// auto-generated by op2.py
//

//user function
__device__ void diff_bflux_gpu( const int *bedgeNum, const double *sJ, const double *nx,
                       const double *ny, const double *s, const double *vis,
                       const double *sigX, const double *sigY, double *flux) {
  int exInd = *bedgeNum * 6;

  for(int i = 0; i < 6; i++) {
    flux[exInd + i] += gaussW_g_cuda[i] * sJ[exInd + i] * 0.5 * (*vis) * s[exInd + i];


  }

}

// CUDA kernel function
__global__ void op_cuda_diff_bflux(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  double *__restrict ind_arg6,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg8_l[18];
    for ( int d=0; d<18; d++ ){
      arg8_l[d] = ZERO_double;
    }
    int map1idx;
    map1idx = opDat1Map[n + set_size * 0];

    //user-supplied kernel call
    diff_bflux_gpu(arg0+n*1,
               ind_arg0+map1idx*18,
               ind_arg1+map1idx*18,
               ind_arg2+map1idx*18,
               ind_arg3+map1idx*18,
               ind_arg4+map1idx*1,
               ind_arg5+map1idx*18,
               ind_arg5+map1idx*18,
               arg8_l);
    atomicAdd(&ind_arg6[0+map1idx*18],arg8_l[0]);
    atomicAdd(&ind_arg6[1+map1idx*18],arg8_l[1]);
    atomicAdd(&ind_arg6[2+map1idx*18],arg8_l[2]);
    atomicAdd(&ind_arg6[3+map1idx*18],arg8_l[3]);
    atomicAdd(&ind_arg6[4+map1idx*18],arg8_l[4]);
    atomicAdd(&ind_arg6[5+map1idx*18],arg8_l[5]);
    atomicAdd(&ind_arg6[6+map1idx*18],arg8_l[6]);
    atomicAdd(&ind_arg6[7+map1idx*18],arg8_l[7]);
    atomicAdd(&ind_arg6[8+map1idx*18],arg8_l[8]);
    atomicAdd(&ind_arg6[9+map1idx*18],arg8_l[9]);
    atomicAdd(&ind_arg6[10+map1idx*18],arg8_l[10]);
    atomicAdd(&ind_arg6[11+map1idx*18],arg8_l[11]);
    atomicAdd(&ind_arg6[12+map1idx*18],arg8_l[12]);
    atomicAdd(&ind_arg6[13+map1idx*18],arg8_l[13]);
    atomicAdd(&ind_arg6[14+map1idx*18],arg8_l[14]);
    atomicAdd(&ind_arg6[15+map1idx*18],arg8_l[15]);
    atomicAdd(&ind_arg6[16+map1idx*18],arg8_l[16]);
    atomicAdd(&ind_arg6[17+map1idx*18],arg8_l[17]);
  }
}


//host stub function
void op_par_loop_diff_bflux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8){

  int nargs = 9;
  op_arg args[9];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(62);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[62].name      = name;
  OP_kernels[62].count    += 1;


  int    ninds   = 7;
  int    inds[9] = {-1,0,1,2,3,4,5,5,6};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: diff_bflux\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_62
      int nthread = OP_BLOCK_SIZE_62;
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
        op_cuda_diff_bflux<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg2.data_d,
        (double *)arg3.data_d,
        (double *)arg4.data_d,
        (double *)arg5.data_d,
        (double *)arg6.data_d,
        (double *)arg8.data_d,
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
  OP_kernels[62].time     += wall_t2 - wall_t1;
}
