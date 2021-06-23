//
// auto-generated by op2.py
//

//user function
__device__ void poisson_mf2_bc_gpu( const double *tol, const int *bedge_type, const int *bedgeNum,
                           const int *d0, const int *d1, const int *d2,
                           const double *mD0, const double *mD1, const double *mD2,
                           const double *sJ, const double *tau, double *op) {
  if(*bedge_type == *d0 || *bedge_type == *d1 || *bedge_type == *d2) {
    for(int j = 0; j < 7 * 15; j++) {
      int indT = (j % 7) * 15 + (j / 7);
      int col  = j % 7;
      int row  = j / 7;
      double val = gaussW_g_cuda[j % 7] * sJ[*bedgeNum * 7 + (j % 7)] * tau[*bedgeNum];
      double mD;
      if(*bedgeNum == 0) {
        val *= gFInterp0_g_cuda[indT];
        mD = mD0[indT];
      } else if(*bedgeNum == 1) {
        val *= gFInterp1_g_cuda[indT];
        mD = mD1[indT];
      } else {
        val *= gFInterp2_g_cuda[indT];
        mD = mD2[indT];
      }
      val -= mD * gaussW_g_cuda[j % 7] * sJ[*bedgeNum * 7 + (j % 7)];
      if(fabs(val) > *tol)
        op[row * 7 + col] += val;
    }
  } else {
    for(int j = 0; j < 7 * 15; j++) {
      int indT = (j % 7) * 15 + (j / 7);
      int col  = j % 7;
      int row  = j / 7;
      double val = gaussW_g_cuda[j % 7] * sJ[*bedgeNum * 7 + (j % 7)];
      if(*bedgeNum == 0) {
        val *= gFInterp0_g_cuda[indT];
      } else if(*bedgeNum == 1) {
        val *= gFInterp1_g_cuda[indT];
      } else {
        val *= gFInterp2_g_cuda[indT];
      }
      if(fabs(val) > *tol)
        op[row * 7 + col] += val;
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_mf2_bc(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const int *__restrict opDat6Map,
  const double *arg0,
  const int *__restrict arg1,
  const int *__restrict arg2,
  const int *arg3,
  const int *arg4,
  const int *arg5,
  double *arg11,
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
    poisson_mf2_bc_gpu(arg0,
                   arg1+n*1,
                   arg2+n*1,
                   arg3,
                   arg4,
                   arg5,
                   ind_arg0+map6idx*105,
                   ind_arg1+map6idx*105,
                   ind_arg2+map6idx*105,
                   ind_arg3+map6idx*21,
                   ind_arg4+map6idx*3,
                   arg11+n*105);
  }
}


//host stub function
void op_par_loop_poisson_mf2_bc(char const *name, op_set set,
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
  op_arg arg11){

  double*arg0h = (double *)arg0.data;
  int*arg3h = (int *)arg3.data;
  int*arg4h = (int *)arg4.data;
  int*arg5h = (int *)arg5.data;
  int nargs = 12;
  op_arg args[12];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(17);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[17].name      = name;
  OP_kernels[17].count    += 1;


  int    ninds   = 5;
  int    inds[12] = {-1,-1,-1,-1,-1,-1,0,1,2,3,4,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_mf2_bc\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(int));
    consts_bytes += ROUND_UP(1*sizeof(int));
    consts_bytes += ROUND_UP(1*sizeof(int));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg0.data   = OP_consts_h + consts_bytes;
    arg0.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg0.data)[d] = arg0h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
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
    arg5.data   = OP_consts_h + consts_bytes;
    arg5.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg5.data)[d] = arg5h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_17
      int nthread = OP_BLOCK_SIZE_17;
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
        op_cuda_poisson_mf2_bc<<<nblocks,nthread>>>(
        (double *)arg6.data_d,
        (double *)arg7.data_d,
        (double *)arg8.data_d,
        (double *)arg9.data_d,
        (double *)arg10.data_d,
        arg6.map_data_d,
        (double*)arg0.data_d,
        (int*)arg1.data_d,
        (int*)arg2.data_d,
        (int*)arg3.data_d,
        (int*)arg4.data_d,
        (int*)arg5.data_d,
        (double*)arg11.data_d,
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
