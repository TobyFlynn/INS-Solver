//
// auto-generated by op2.py
//

//user function
__device__ void viscosity_bc_gpu( const int *bedge_type, const int *bedgeNum,
                         const double *t, const int *problem, const double *x,
                         const double *y, const double *nx, const double *ny,
                         const double *nu, double *exQ0, double *exQ1) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 7;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 7;
  }

  const double PI = 3.141592653589793238463;

  if(*problem == 0) {
    if(*bedge_type == 0) {

      for(int i = 0; i < 7; i++) {
        double y1 = y[exInd + i];
        exQ0[exInd + i] += pow(1.0, -2.0) * sin((PI * (*t)) / 8.0) * 6.0 * y1 * (1.0 - y1);
      }
    } else if(*bedge_type == 1) {






    } else {

    }
  } else {
    if(*bedge_type == 0) {

      for(int i = 0; i < 7; i++) {
        double y1 = y[exInd + i];
        double x1 = x[exInd + i];
        exQ0[exInd + i] += -sin(2.0 * PI * y1) * exp(-nu[exInd + i] * 4.0 * PI * PI * *t);
        exQ1[exInd + i] += sin(2.0 * PI * x1) * exp(-nu[exInd + i] * 4.0 * PI * PI * *t);
      }
    }

    if(*bedge_type == 1) {

      for(int i = 0; i < 7; i++) {
        double y1  = y[exInd + i];
        double x1  = x[exInd + i];
        double ny1 = ny[exInd + i];
        double nx1 = nx[exInd + i];
        exQ0[exInd + i] += ny1 * 2.0 * PI * (-cos(2.0 * PI * y1)) * exp(-nu[exInd + i] * 4.0 * PI * PI * *t);
        exQ1[exInd + i] += nx1 * 2.0 * PI * cos(2.0 * PI * x1) * exp(-nu[exInd + i] * 4.0 * PI * PI * *t);
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_viscosity_bc(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  double *__restrict ind_arg5,
  double *__restrict ind_arg6,
  const int *__restrict opDat4Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  const double *arg2,
  const int *arg3,
  int start,
  int end,
  int   set_size) {
  double arg9_l[21];
  double arg10_l[21];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg9_l[21];
    for ( int d=0; d<21; d++ ){
      arg9_l[d] = ZERO_double;
    }
    double arg10_l[21];
    for ( int d=0; d<21; d++ ){
      arg10_l[d] = ZERO_double;
    }
    int map4idx;
    map4idx = opDat4Map[n + set_size * 0];

    //user-supplied kernel call
    viscosity_bc_gpu(arg0+n*1,
                 arg1+n*1,
                 arg2,
                 arg3,
                 ind_arg0+map4idx*21,
                 ind_arg1+map4idx*21,
                 ind_arg2+map4idx*21,
                 ind_arg3+map4idx*21,
                 ind_arg4+map4idx*21,
                 arg9_l,
                 arg10_l);
    atomicAdd(&ind_arg5[0+map4idx*21],arg9_l[0]);
    atomicAdd(&ind_arg5[1+map4idx*21],arg9_l[1]);
    atomicAdd(&ind_arg5[2+map4idx*21],arg9_l[2]);
    atomicAdd(&ind_arg5[3+map4idx*21],arg9_l[3]);
    atomicAdd(&ind_arg5[4+map4idx*21],arg9_l[4]);
    atomicAdd(&ind_arg5[5+map4idx*21],arg9_l[5]);
    atomicAdd(&ind_arg5[6+map4idx*21],arg9_l[6]);
    atomicAdd(&ind_arg5[7+map4idx*21],arg9_l[7]);
    atomicAdd(&ind_arg5[8+map4idx*21],arg9_l[8]);
    atomicAdd(&ind_arg5[9+map4idx*21],arg9_l[9]);
    atomicAdd(&ind_arg5[10+map4idx*21],arg9_l[10]);
    atomicAdd(&ind_arg5[11+map4idx*21],arg9_l[11]);
    atomicAdd(&ind_arg5[12+map4idx*21],arg9_l[12]);
    atomicAdd(&ind_arg5[13+map4idx*21],arg9_l[13]);
    atomicAdd(&ind_arg5[14+map4idx*21],arg9_l[14]);
    atomicAdd(&ind_arg5[15+map4idx*21],arg9_l[15]);
    atomicAdd(&ind_arg5[16+map4idx*21],arg9_l[16]);
    atomicAdd(&ind_arg5[17+map4idx*21],arg9_l[17]);
    atomicAdd(&ind_arg5[18+map4idx*21],arg9_l[18]);
    atomicAdd(&ind_arg5[19+map4idx*21],arg9_l[19]);
    atomicAdd(&ind_arg5[20+map4idx*21],arg9_l[20]);
    atomicAdd(&ind_arg6[0+map4idx*21],arg10_l[0]);
    atomicAdd(&ind_arg6[1+map4idx*21],arg10_l[1]);
    atomicAdd(&ind_arg6[2+map4idx*21],arg10_l[2]);
    atomicAdd(&ind_arg6[3+map4idx*21],arg10_l[3]);
    atomicAdd(&ind_arg6[4+map4idx*21],arg10_l[4]);
    atomicAdd(&ind_arg6[5+map4idx*21],arg10_l[5]);
    atomicAdd(&ind_arg6[6+map4idx*21],arg10_l[6]);
    atomicAdd(&ind_arg6[7+map4idx*21],arg10_l[7]);
    atomicAdd(&ind_arg6[8+map4idx*21],arg10_l[8]);
    atomicAdd(&ind_arg6[9+map4idx*21],arg10_l[9]);
    atomicAdd(&ind_arg6[10+map4idx*21],arg10_l[10]);
    atomicAdd(&ind_arg6[11+map4idx*21],arg10_l[11]);
    atomicAdd(&ind_arg6[12+map4idx*21],arg10_l[12]);
    atomicAdd(&ind_arg6[13+map4idx*21],arg10_l[13]);
    atomicAdd(&ind_arg6[14+map4idx*21],arg10_l[14]);
    atomicAdd(&ind_arg6[15+map4idx*21],arg10_l[15]);
    atomicAdd(&ind_arg6[16+map4idx*21],arg10_l[16]);
    atomicAdd(&ind_arg6[17+map4idx*21],arg10_l[17]);
    atomicAdd(&ind_arg6[18+map4idx*21],arg10_l[18]);
    atomicAdd(&ind_arg6[19+map4idx*21],arg10_l[19]);
    atomicAdd(&ind_arg6[20+map4idx*21],arg10_l[20]);
  }
}


//host stub function
void op_par_loop_viscosity_bc(char const *name, op_set set,
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

  double*arg2h = (double *)arg2.data;
  int*arg3h = (int *)arg3.data;
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
  op_timing_realloc(36);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[36].name      = name;
  OP_kernels[36].count    += 1;


  int    ninds   = 7;
  int    inds[11] = {-1,-1,-1,-1,0,1,2,3,4,5,6};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: viscosity_bc\n");
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
    #ifdef OP_BLOCK_SIZE_36
      int nthread = OP_BLOCK_SIZE_36;
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
        op_cuda_viscosity_bc<<<nblocks,nthread>>>(
        (double *)arg4.data_d,
        (double *)arg5.data_d,
        (double *)arg6.data_d,
        (double *)arg7.data_d,
        (double *)arg8.data_d,
        (double *)arg9.data_d,
        (double *)arg10.data_d,
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
  OP_kernels[36].time     += wall_t2 - wall_t1;
}
