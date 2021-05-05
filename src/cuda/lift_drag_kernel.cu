//
// auto-generated by op2.py
//

//user function
__device__ void lift_drag_gpu( const int *bedge_type, const int *bedgeNum, const double *p,
                      const double *dQ0dx, const double *dQ0dy,
                      const double *dQ1dx, const double *dQ1dy,
                      const double *nx, const double *ny, const double *sJ,
                      double *cd, double *cl) {
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

  if(*bedge_type == 2) {
    for(int i = 0; i < 5; i++) {
      *cd += lift_drag_vec_cuda[i] * sJ[exInd + i] * (-p[fmask[i]] * nx[exInd + i] + nu_cuda * (nx[exInd + i] * 2.0 * dQ0dx[fmask[i]] + ny[exInd + i] * (dQ1dx[fmask[i]] + dQ0dy[fmask[i]])));
      *cl += lift_drag_vec_cuda[i] * sJ[exInd + i] * (-p[fmask[i]] * ny[exInd + i] + nu_cuda * (nx[exInd + i] * (dQ1dx[fmask[i]] + dQ0dy[fmask[i]]) + ny[exInd + i] * 2.0 * dQ1dy[fmask[i]]));
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_lift_drag(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  const double *__restrict ind_arg6,
  const double *__restrict ind_arg7,
  const int *__restrict opDat2Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  double *arg10,
  double *arg11,
  int start,
  int end,
  int   set_size) {
  double arg10_l[1];
  for ( int d=0; d<1; d++ ){
    arg10_l[d]=ZERO_double;
  }
  double arg11_l[1];
  for ( int d=0; d<1; d++ ){
    arg11_l[d]=ZERO_double;
  }
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    int map2idx;
    map2idx = opDat2Map[n + set_size * 0];

    //user-supplied kernel call
    lift_drag_gpu(arg0+n*1,
              arg1+n*1,
              ind_arg0+map2idx*15,
              ind_arg1+map2idx*15,
              ind_arg2+map2idx*15,
              ind_arg3+map2idx*15,
              ind_arg4+map2idx*15,
              ind_arg5+map2idx*15,
              ind_arg6+map2idx*15,
              ind_arg7+map2idx*15,
              arg10_l,
              arg11_l);
  }

  //global reductions

  for ( int d=0; d<1; d++ ){
    op_reduction<OP_INC>(&arg10[d+blockIdx.x*1],arg10_l[d]);
  }
  for ( int d=0; d<1; d++ ){
    op_reduction<OP_INC>(&arg11[d+blockIdx.x*1],arg11_l[d]);
  }
}


//host stub function
void op_par_loop_lift_drag(char const *name, op_set set,
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

  double*arg10h = (double *)arg10.data;
  double*arg11h = (double *)arg11.data;
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
  op_timing_realloc(57);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[57].name      = name;
  OP_kernels[57].count    += 1;


  int    ninds   = 8;
  int    inds[12] = {-1,-1,0,1,2,3,4,5,6,7,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: lift_drag\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_57
      int nthread = OP_BLOCK_SIZE_57;
    #else
      int nthread = OP_block_size;
    #endif

    //transfer global reduction data to GPU
    int maxblocks = (MAX(set->core_size, set->size+set->exec_size-set->core_size)-1)/nthread+1;
    int reduct_bytes = 0;
    int reduct_size  = 0;
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    reduct_size   = MAX(reduct_size,sizeof(double));
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    reduct_size   = MAX(reduct_size,sizeof(double));
    reallocReductArrays(reduct_bytes);
    reduct_bytes = 0;
    arg10.data   = OP_reduct_h + reduct_bytes;
    arg10.data_d = OP_reduct_d + reduct_bytes;
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        ((double *)arg10.data)[d+b*1] = ZERO_double;
      }
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    arg11.data   = OP_reduct_h + reduct_bytes;
    arg11.data_d = OP_reduct_d + reduct_bytes;
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        ((double *)arg11.data)[d+b*1] = ZERO_double;
      }
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    mvReductArraysToDevice(reduct_bytes);

    for ( int round=0; round<3; round++ ){
      if (round==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = round==0 ? 0 : (round==1 ? set->core_size : set->size);
      int end = round==0 ? set->core_size : (round==1? set->size :  set->size + set->exec_size);
      if (end-start>0) {
        int nblocks = (end-start-1)/nthread+1;
        int nshared = reduct_size*nthread;
        op_cuda_lift_drag<<<nblocks,nthread,nshared>>>(
        (double *)arg2.data_d,
        (double *)arg3.data_d,
        (double *)arg4.data_d,
        (double *)arg5.data_d,
        (double *)arg6.data_d,
        (double *)arg7.data_d,
        (double *)arg8.data_d,
        (double *)arg9.data_d,
        arg2.map_data_d,
        (int*)arg0.data_d,
        (int*)arg1.data_d,
        (double*)arg10.data_d,
        (double*)arg11.data_d,
        start,end,set->size+set->exec_size);
      }
      if (round==1) mvReductArraysToHost(reduct_bytes);
    }
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        arg10h[d] = arg10h[d] + ((double *)arg10.data)[d+b*1];
      }
    }
    arg10.data = (char *)arg10h;
    op_mpi_reduce(&arg10,arg10h);
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        arg11h[d] = arg11h[d] + ((double *)arg11.data)[d+b*1];
      }
    }
    arg11.data = (char *)arg11h;
    op_mpi_reduce(&arg11,arg11h);
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[57].time     += wall_t2 - wall_t1;
}
