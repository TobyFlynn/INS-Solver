//
// auto-generated by op2.py
//

//user function
__device__ void init_gauss_grad5_gpu( const int *edgeNum,
                             const double *nxL, const double *nxR,
                             const double *nyL, const double *nyR,
                             const double *Dx0L, const double *Dx0R,
                             const double *Dy0L, const double *Dy0R,
                             const double *Dx1L, const double *Dx1R,
                             const double *Dy1L, const double *Dy1R,
                             const double *Dx2L, const double *Dx2R,
                             const double *Dy2L, const double *Dy2R,
                             double *dL, double *dR) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];

  const double *DxL, *DyL, *DxR, *DyR;

  if(edgeL == 0) {
    DxR = Dx0L;
    DyR = Dy0L;
  } else if(edgeL == 1) {
    DxR = Dx1L;
    DyR = Dy1L;
  } else {
    DxR = Dx2L;
    DyR = Dy2L;
  }

  if(edgeR == 0) {
    DxL = Dx0R;
    DyL = Dy0R;
  } else if(edgeR == 1) {
    DxL = Dx1R;
    DyL = Dy1R;
  } else {
    DxL = Dx2R;
    DyL = Dy2R;
  }

  for(int j = 0; j < 10; j++) {
    for(int i = 0; i < 6; i++) {
      int ind  = j * 6 + i;
      int indL = edgeL * 6 + i;
      int indR = edgeR * 6 + i;

      dL[ind] = nxL[indL] * DxL[ind] + nyL[indL] * DyL[ind];
      dR[ind] = nxR[indR] * DxR[ind] + nyR[indR] * DyR[ind];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_init_gauss_grad5(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  const double *__restrict ind_arg6,
  const double *__restrict ind_arg7,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  double *arg17,
  double *arg18,
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

    //user-supplied kernel call
    init_gauss_grad5_gpu(arg0+n*2,
                     ind_arg0+map1idx*18,
                     ind_arg0+map2idx*18,
                     ind_arg1+map1idx*18,
                     ind_arg1+map2idx*18,
                     ind_arg2+map1idx*60,
                     ind_arg2+map2idx*60,
                     ind_arg3+map1idx*60,
                     ind_arg3+map2idx*60,
                     ind_arg4+map1idx*60,
                     ind_arg4+map2idx*60,
                     ind_arg5+map1idx*60,
                     ind_arg5+map2idx*60,
                     ind_arg6+map1idx*60,
                     ind_arg6+map2idx*60,
                     ind_arg7+map1idx*60,
                     ind_arg7+map2idx*60,
                     arg17+n*60,
                     arg18+n*60);
  }
}


//host stub function
void op_par_loop_init_gauss_grad5(char const *name, op_set set,
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
  op_arg arg14,
  op_arg arg15,
  op_arg arg16,
  op_arg arg17,
  op_arg arg18){

  int nargs = 19;
  op_arg args[19];

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
  args[15] = arg15;
  args[16] = arg16;
  args[17] = arg17;
  args[18] = arg18;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(7);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[7].name      = name;
  OP_kernels[7].count    += 1;


  int    ninds   = 8;
  int    inds[19] = {-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: init_gauss_grad5\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_7
      int nthread = OP_BLOCK_SIZE_7;
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
        op_cuda_init_gauss_grad5<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg3.data_d,
        (double *)arg5.data_d,
        (double *)arg7.data_d,
        (double *)arg9.data_d,
        (double *)arg11.data_d,
        (double *)arg13.data_d,
        (double *)arg15.data_d,
        arg1.map_data_d,
        (int*)arg0.data_d,
        (double*)arg17.data_d,
        (double*)arg18.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[7].time     += wall_t2 - wall_t1;
}
