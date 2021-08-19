//
// auto-generated by op2.py
//

//user function
__device__ void poisson_op5_gpu( const int *edgeType, const int *edgeNum,
                        const int *d0, const int *d1, const int *d2,
                        const double *mD, const double *sJ, const double *h,
                        const double *gFactor, const double *factor,
                        double *op) {

  const double *gVM;
  if(*edgeNum == 0) {
    gVM = gFInterp0_g_cuda;
  } else if(*edgeNum == 1) {
    gVM = gFInterp1_g_cuda;
  } else {
    gVM = gFInterp2_g_cuda;
  }

  for(int i = 0; i < 6 * 10; i++) {
    op[i] = 0.0;
  }

  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2) {


    for(int i = 0; i < 6 * 10; i++) {
      int indT = (i % 6) * 10 + i / 6;
      int indSJ = *edgeNum * 6 + (i % 6);
      op[i] = gVM[indT] * gaussW_g_cuda[i % 6] * sJ[indSJ];
    }
  } else {

    double tauA[6];
    double maxTau = 0.0;
    for(int i = 0; i < 6; i++) {
      int ind = *edgeNum  * 6 + i;

      tauA[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * (*h * gFactor[ind]);


    }





    for(int i = 0; i < 6 * 10; i++) {
      int indT = (i % 6) * 10 + i / 6;
      int indSJ = *edgeNum * 6 + (i % 6);
      int indFactor = (i / 6);

      op[i] = gVM[indT] * gaussW_g_cuda[i % 6] * sJ[indSJ] * tauA[i % 6]
              - factor[indFactor] * mD[indT] * gaussW_g_cuda[i % 6] * sJ[indSJ];


    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_op5(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const int *__restrict opDat6Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  const int *arg2,
  const int *arg3,
  const int *arg4,
  const double *__restrict arg5,
  double *arg10,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    int map6idx;
    map6idx = opDat6Map[n + set_size * 0];

    //user-supplied kernel call
    poisson_op5_gpu(arg0+n*1,
                arg1+n*1,
                arg2,
                arg3,
                arg4,
                arg5+n*60,
                ind_arg0+map6idx*18,
                ind_arg1+map6idx*1,
                ind_arg2+map6idx*18,
                ind_arg3+map6idx*10,
                arg10+n*60);
  }
}


//host stub function
void op_par_loop_poisson_op5(char const *name, op_set set,
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
  op_arg arg10){

  int*arg2h = (int *)arg2.data;
  int*arg3h = (int *)arg3.data;
  int*arg4h = (int *)arg4.data;
  int nargs = 11;
  op_arg args[11];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(20);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[20].name      = name;
  OP_kernels[20].count    += 1;


  int    ninds   = 4;
  int    inds[11] = {-1,-1,-1,-1,-1,-1,0,1,2,3,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_op5\n");
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
    #ifdef OP_BLOCK_SIZE_20
      int nthread = OP_BLOCK_SIZE_20;
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
        op_cuda_poisson_op5<<<nblocks,nthread>>>(
        (double *)arg6.data_d,
        (double *)arg7.data_d,
        (double *)arg8.data_d,
        (double *)arg9.data_d,
        arg6.map_data_d,
        (int*)arg0.data_d,
        (int*)arg1.data_d,
        (int*)arg2.data_d,
        (int*)arg3.data_d,
        (int*)arg4.data_d,
        (double*)arg5.data_d,
        (double*)arg10.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[20].time     += wall_t2 - wall_t1;
}
