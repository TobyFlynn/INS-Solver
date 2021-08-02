//
// auto-generated by op2.py
//

//user function
__device__ void pressure_bc_gpu( const int *bedge_type, const int *bedgeNum,
                        const double *t, const int *problem, const double *x,
                        const double *y, const double *nx, const double *ny,
                        const double *nu, const double *rho, const double *N0, const double *N1,
                        const double *gradCurlVel0, const double *gradCurlVel1,
                        double *dPdN) {
  int exInd = *bedgeNum * 4;
  int *fmask = &FMASK_cuda[*bedgeNum * 4];

  const double PI = 3.141592653589793238463;

  if(*problem == 0) {
    if(*bedge_type == 0 || *bedge_type == 2 || *bedge_type == 3) {

      for(int i = 0; i < 4; i++) {
        int fInd = fmask[i];


        double res1 = -N0[fInd] - gradCurlVel1[fInd] / (reynolds_cuda * rho[fInd]);
        double res2 = -N1[fInd] + gradCurlVel0[fInd] / (reynolds_cuda * rho[fInd]);
        dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;
      }
    }

    if(*bedge_type == 0) {

      for(int i = 0; i < 4; i++) {
        double y1 = y[fmask[i]];
        double bcdUndt = -pow(1.0, -2.0) * (PI/8.0) * cos((PI * *t) / 8.0) * 6.0 * y1 * (1.0 - y1);
        dPdN[exInd + i] -= bcdUndt;
      }
    }
  } else {
    if(*bedge_type == 0) {

      for(int i = 0; i < 4; i++) {
        int fInd = fmask[i];
        double res1 = -N0[fInd] - nu[fInd] * gradCurlVel1[fInd];
        double res2 = -N1[fInd] + nu[fInd] * gradCurlVel0[fInd];
        dPdN[exInd + i] += nx[exInd + i] * res1 + ny[exInd + i] * res2;

        double y1 = y[fmask[i]];
        double x1 = x[fmask[i]];
        double nx1 = nx[exInd + i];
        double ny1 = ny[exInd + i];
        double bcdUndt = -nu[fInd] * 4.0 * PI * PI * (-nx1 * sin(2.0 * PI * y1) + ny1 * sin(2.0 * PI * x1))
                          * exp(-nu[fInd] * 4.0 * PI * PI * *t);
        dPdN[exInd + i] -= bcdUndt;
      }
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
  const double *__restrict ind_arg6,
  const double *__restrict ind_arg7,
  const double *__restrict ind_arg8,
  const double *__restrict ind_arg9,
  double *__restrict ind_arg10,
  const int *__restrict opDat4Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  const double *arg2,
  const int *arg3,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg14_l[12];
    for ( int d=0; d<12; d++ ){
      arg14_l[d] = ZERO_double;
    }
    int map4idx;
    map4idx = opDat4Map[n + set_size * 0];

    //user-supplied kernel call
    pressure_bc_gpu(arg0+n*1,
                arg1+n*1,
                arg2,
                arg3,
                ind_arg0+map4idx*10,
                ind_arg1+map4idx*10,
                ind_arg2+map4idx*12,
                ind_arg3+map4idx*12,
                ind_arg4+map4idx*10,
                ind_arg5+map4idx*10,
                ind_arg6+map4idx*10,
                ind_arg7+map4idx*10,
                ind_arg8+map4idx*10,
                ind_arg9+map4idx*10,
                arg14_l);
    atomicAdd(&ind_arg10[0+map4idx*12],arg14_l[0]);
    atomicAdd(&ind_arg10[1+map4idx*12],arg14_l[1]);
    atomicAdd(&ind_arg10[2+map4idx*12],arg14_l[2]);
    atomicAdd(&ind_arg10[3+map4idx*12],arg14_l[3]);
    atomicAdd(&ind_arg10[4+map4idx*12],arg14_l[4]);
    atomicAdd(&ind_arg10[5+map4idx*12],arg14_l[5]);
    atomicAdd(&ind_arg10[6+map4idx*12],arg14_l[6]);
    atomicAdd(&ind_arg10[7+map4idx*12],arg14_l[7]);
    atomicAdd(&ind_arg10[8+map4idx*12],arg14_l[8]);
    atomicAdd(&ind_arg10[9+map4idx*12],arg14_l[9]);
    atomicAdd(&ind_arg10[10+map4idx*12],arg14_l[10]);
    atomicAdd(&ind_arg10[11+map4idx*12],arg14_l[11]);
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
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14){

  double*arg2h = (double *)arg2.data;
  int*arg3h = (int *)arg3.data;
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
  op_timing_realloc(35);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[35].name      = name;
  OP_kernels[35].count    += 1;


  int    ninds   = 11;
  int    inds[15] = {-1,-1,-1,-1,0,1,2,3,4,5,6,7,8,9,10};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_bc\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
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
    #ifdef OP_BLOCK_SIZE_35
      int nthread = OP_BLOCK_SIZE_35;
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
        op_cuda_pressure_bc<<<nblocks,nthread>>>(
        (double *)arg4.data_d,
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
  OP_kernels[35].time     += wall_t2 - wall_t1;
}
