//
// auto-generated by op2.py
//

//user function
__device__ void pressure_bc_gpu( const int *bedge_type, const int *bedgeNum,
                        const double *nx, const double *ny,
                        const double *N0, const double *N1,
                        const double *gradCurlVel0, const double *gradCurlVel1,
                        double *dPdN) {
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

  if(*bedge_type == 0 || *bedge_type == 2) {

    for(int i = 0; i < 5; i++) {
      int fInd = fmask[i];
      double res1 = -N0[fInd] - nu_cuda * gradCurlVel1[fInd];
      double res2 = -N1[fInd] + nu_cuda * gradCurlVel0[fInd];
      dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;
    }
  }

  if(*bedge_type == 0) {


    double bcdUndt = -1.0;
    for(int i = 0; i < 5; i++) {
      dPdN[exInd + i] -= bcdUndt;
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_pressure_bc(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  double *__restrict ind_arg6,
  const int *__restrict opDat2Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  int start,
  int end,
  int   set_size) {
  double arg8_l[15];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg8_l[15];
    for ( int d=0; d<15; d++ ){
      arg8_l[d] = ZERO_double;
    }
    int map2idx;
    map2idx = opDat2Map[n + set_size * 0];

    //user-supplied kernel call
    pressure_bc_gpu(arg0+n*1,
                arg1+n*1,
                ind_arg0+map2idx*15,
                ind_arg1+map2idx*15,
                ind_arg2+map2idx*15,
                ind_arg3+map2idx*15,
                ind_arg4+map2idx*15,
                ind_arg5+map2idx*15,
                arg8_l);
    atomicAdd(&ind_arg6[0+map2idx*15],arg8_l[0]);
    atomicAdd(&ind_arg6[1+map2idx*15],arg8_l[1]);
    atomicAdd(&ind_arg6[2+map2idx*15],arg8_l[2]);
    atomicAdd(&ind_arg6[3+map2idx*15],arg8_l[3]);
    atomicAdd(&ind_arg6[4+map2idx*15],arg8_l[4]);
    atomicAdd(&ind_arg6[5+map2idx*15],arg8_l[5]);
    atomicAdd(&ind_arg6[6+map2idx*15],arg8_l[6]);
    atomicAdd(&ind_arg6[7+map2idx*15],arg8_l[7]);
    atomicAdd(&ind_arg6[8+map2idx*15],arg8_l[8]);
    atomicAdd(&ind_arg6[9+map2idx*15],arg8_l[9]);
    atomicAdd(&ind_arg6[10+map2idx*15],arg8_l[10]);
    atomicAdd(&ind_arg6[11+map2idx*15],arg8_l[11]);
    atomicAdd(&ind_arg6[12+map2idx*15],arg8_l[12]);
    atomicAdd(&ind_arg6[13+map2idx*15],arg8_l[13]);
    atomicAdd(&ind_arg6[14+map2idx*15],arg8_l[14]);
  }
}


//host stub function
void op_par_loop_pressure_bc(char const *name, op_set set,
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
  op_timing_realloc(7);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[7].name      = name;
  OP_kernels[7].count    += 1;


  int    ninds   = 7;
  int    inds[9] = {-1,-1,0,1,2,3,4,5,6};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_bc\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_7
      int nthread = OP_BLOCK_SIZE_7;
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
        op_cuda_pressure_bc<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg3.data_d,
        (double *)arg4.data_d,
        (double *)arg5.data_d,
        (double *)arg6.data_d,
        (double *)arg7.data_d,
        (double *)arg8.data_d,
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
  OP_kernels[7].time     += wall_t2 - wall_t1;
}
