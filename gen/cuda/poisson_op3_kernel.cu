//
// auto-generated by op2.py
//

//user function
__device__ void poisson_op3_gpu( const int *edgeType, const int *edgeNum,
                        const int *d0, const int *d1, const int *d2,
                        const double *mD0, const double *mD1,
                        const double *mD2, const double *sJ,
                        const double *h, const double *gFactor,
                        const double *factor, double *op1) {
  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2)
    return;


  const double *mD, *gVM;
  if(*edgeNum == 0) {
    mD  = mD0;
    gVM = gFInterp0_g_cuda;
  } else if(*edgeNum == 1) {
    mD  = mD1;
    gVM = gFInterp1_g_cuda;
  } else {
    mD  = mD2;
    gVM = gFInterp2_g_cuda;
  }


  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      int c_ind = i * 3 + j;
      for(int k = 0; k < 3; k++) {

        int b_ind = k * 3 + j;

        int ind = i * 3 + k;
        int a_ind = ((ind * 3) % (3 * 3)) + (ind / 3);

        int factors_ind = *edgeNum * 3 + k;

        op1[c_ind] += -0.5 * gVM[a_ind] * gaussW_g_cuda[k] * sJ[factors_ind]
                      * gFactor[factors_ind] * mD[b_ind];
      }
    }
  }


  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      int c_ind = i * 3 + j;
      for(int k = 0; k < 3; k++) {

        int b_ind = k * 3 + j;

        int ind = i * 3 + k;
        int a_ind = ((ind * 3) % (3 * 3)) + (ind / 3);

        int factors_ind = *edgeNum * 3 + k;

        op1[c_ind] += -factor[i] * mD[a_ind] * gaussW_g_cuda[k]
                      * sJ[factors_ind] * gVM[b_ind];
      }
    }
  }


  double tauA[3];
  for(int i = 0; i < 3; i++) {
    int ind = *edgeNum  * 3 + i;
    tauA[i] = 100 * 0.5 * 5 * 6 * (*h * gFactor[ind]);

  }



  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      int c_ind = i * 3 + j;
      for(int k = 0; k < 3; k++) {

        int b_ind = k * 3 + j;

        int ind = i * 3 + k;
        int a_ind = ((ind * 3) % (3 * 3)) + (ind / 3);

        int factors_ind = *edgeNum * 3 + k;

        op1[c_ind] += gVM[a_ind] * gaussW_g_cuda[k] * sJ[factors_ind]
                      * tauA[k] * gVM[b_ind];
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_op3(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  const double *__restrict ind_arg6,
  double *__restrict ind_arg7,
  const int *__restrict opDat5Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  const int *arg2,
  const int *arg3,
  const int *arg4,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg12_l[9];
    for ( int d=0; d<9; d++ ){
      arg12_l[d] = ZERO_double;
    }
    int map5idx;
    map5idx = opDat5Map[n + set_size * 0];

    //user-supplied kernel call
    poisson_op3_gpu(arg0+n*1,
                arg1+n*1,
                arg2,
                arg3,
                arg4,
                ind_arg0+map5idx*9,
                ind_arg1+map5idx*9,
                ind_arg2+map5idx*9,
                ind_arg3+map5idx*9,
                ind_arg4+map5idx*1,
                ind_arg5+map5idx*9,
                ind_arg6+map5idx*3,
                arg12_l);
    atomicAdd(&ind_arg7[0+map5idx*9],arg12_l[0]);
    atomicAdd(&ind_arg7[1+map5idx*9],arg12_l[1]);
    atomicAdd(&ind_arg7[2+map5idx*9],arg12_l[2]);
    atomicAdd(&ind_arg7[3+map5idx*9],arg12_l[3]);
    atomicAdd(&ind_arg7[4+map5idx*9],arg12_l[4]);
    atomicAdd(&ind_arg7[5+map5idx*9],arg12_l[5]);
    atomicAdd(&ind_arg7[6+map5idx*9],arg12_l[6]);
    atomicAdd(&ind_arg7[7+map5idx*9],arg12_l[7]);
    atomicAdd(&ind_arg7[8+map5idx*9],arg12_l[8]);
  }
}


//host stub function
void op_par_loop_poisson_op3(char const *name, op_set set,
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
  op_arg arg12){

  int*arg2h = (int *)arg2.data;
  int*arg3h = (int *)arg3.data;
  int*arg4h = (int *)arg4.data;
  int nargs = 13;
  op_arg args[13];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(21);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[21].name      = name;
  OP_kernels[21].count    += 1;


  int    ninds   = 8;
  int    inds[13] = {-1,-1,-1,-1,-1,0,1,2,3,4,5,6,7};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_op3\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
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
    #ifdef OP_BLOCK_SIZE_21
      int nthread = OP_BLOCK_SIZE_21;
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
        op_cuda_poisson_op3<<<nblocks,nthread>>>(
        (double *)arg5.data_d,
        (double *)arg6.data_d,
        (double *)arg7.data_d,
        (double *)arg8.data_d,
        (double *)arg9.data_d,
        (double *)arg10.data_d,
        (double *)arg11.data_d,
        (double *)arg12.data_d,
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
  OP_kernels[21].time     += wall_t2 - wall_t1;
}
