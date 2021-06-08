//
// auto-generated by op2.py
//

//user function
__device__ void pressure_bc2_gpu( const int *bedge_type, const int *bedgeNum,
                         const double *t, const int *problem, const double *x,
                         const double *y, double *prBC) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 7;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 7;
  }

  const double PI = 3.141592653589793238463;

  if(*problem == 1) {
    if(*bedge_type == 1) {

      for(int i = 0; i < 7; i++) {
        double y1 = y[exInd + i];
        double x1 = x[exInd + i];
        prBC[exInd + i] += -cos(2.0 * PI * x1) * cos(2.0 * PI * y1) * exp(-nu_cuda * 8.0 * PI * PI * *t);
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_pressure_bc2(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  double *__restrict ind_arg2,
  const int *__restrict opDat4Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  const double *arg2,
  const int *arg3,
  int start,
  int end,
  int   set_size) {
  double arg6_l[21];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg6_l[21];
    for ( int d=0; d<21; d++ ){
      arg6_l[d] = ZERO_double;
    }
    int map4idx;
    map4idx = opDat4Map[n + set_size * 0];

    //user-supplied kernel call
    pressure_bc2_gpu(arg0+n*1,
                 arg1+n*1,
                 arg2,
                 arg3,
                 ind_arg0+map4idx*21,
                 ind_arg1+map4idx*21,
                 arg6_l);
    atomicAdd(&ind_arg2[0+map4idx*21],arg6_l[0]);
    atomicAdd(&ind_arg2[1+map4idx*21],arg6_l[1]);
    atomicAdd(&ind_arg2[2+map4idx*21],arg6_l[2]);
    atomicAdd(&ind_arg2[3+map4idx*21],arg6_l[3]);
    atomicAdd(&ind_arg2[4+map4idx*21],arg6_l[4]);
    atomicAdd(&ind_arg2[5+map4idx*21],arg6_l[5]);
    atomicAdd(&ind_arg2[6+map4idx*21],arg6_l[6]);
    atomicAdd(&ind_arg2[7+map4idx*21],arg6_l[7]);
    atomicAdd(&ind_arg2[8+map4idx*21],arg6_l[8]);
    atomicAdd(&ind_arg2[9+map4idx*21],arg6_l[9]);
    atomicAdd(&ind_arg2[10+map4idx*21],arg6_l[10]);
    atomicAdd(&ind_arg2[11+map4idx*21],arg6_l[11]);
    atomicAdd(&ind_arg2[12+map4idx*21],arg6_l[12]);
    atomicAdd(&ind_arg2[13+map4idx*21],arg6_l[13]);
    atomicAdd(&ind_arg2[14+map4idx*21],arg6_l[14]);
    atomicAdd(&ind_arg2[15+map4idx*21],arg6_l[15]);
    atomicAdd(&ind_arg2[16+map4idx*21],arg6_l[16]);
    atomicAdd(&ind_arg2[17+map4idx*21],arg6_l[17]);
    atomicAdd(&ind_arg2[18+map4idx*21],arg6_l[18]);
    atomicAdd(&ind_arg2[19+map4idx*21],arg6_l[19]);
    atomicAdd(&ind_arg2[20+map4idx*21],arg6_l[20]);
  }
}


//host stub function
void op_par_loop_pressure_bc2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6){

  double*arg2h = (double *)arg2.data;
  int*arg3h = (int *)arg3.data;
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
  op_timing_realloc(40);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[40].name      = name;
  OP_kernels[40].count    += 1;


  int    ninds   = 3;
  int    inds[7] = {-1,-1,-1,-1,0,1,2};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_bc2\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(int));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg2.data   = OP_consts_h + consts_bytes;
    arg2.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg2.data)[d] = arg2h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    arg3.data   = OP_consts_h + consts_bytes;
    arg3.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg3.data)[d] = arg3h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_40
      int nthread = OP_BLOCK_SIZE_40;
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
        op_cuda_pressure_bc2<<<nblocks,nthread>>>(
        (double *)arg4.data_d,
        (double *)arg5.data_d,
        (double *)arg6.data_d,
        arg4.map_data_d,
        (int*)arg0.data_d,
        (int*)arg1.data_d,
        (double*)arg2.data_d,
        (int*)arg3.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[40].time     += wall_t2 - wall_t1;
}
