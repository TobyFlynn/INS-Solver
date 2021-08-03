//
// auto-generated by op2.py
//

//user function
__device__ void pressure_bc2_gpu( const int *bedge_type, const int *bedgeNum,
                         const double *t, const int *problem, const double *x,
                         const double *y, const double *nu, double *prBC) {
  int exInd = *bedgeNum * 6;

  const double PI = 3.141592653589793238463;

  if(*problem == 1) {
    if(*bedge_type == 1) {

      for(int i = 0; i < 6; i++) {
        double y1 = y[exInd + i];
        double x1 = x[exInd + i];
        prBC[exInd + i] += -cos(2.0 * PI * x1) * cos(2.0 * PI * y1) * exp(-nu[exInd + i] * 8.0 * PI * PI * *t);
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_pressure_bc2(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  double *__restrict ind_arg3,
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
    double arg7_l[18];
    for ( int d=0; d<18; d++ ){
      arg7_l[d] = ZERO_double;
    }
    int map4idx;
    map4idx = opDat4Map[n + set_size * 0];

    //user-supplied kernel call
    pressure_bc2_gpu(arg0+n*1,
                 arg1+n*1,
                 arg2,
                 arg3,
                 ind_arg0+map4idx*18,
                 ind_arg1+map4idx*18,
                 ind_arg2+map4idx*18,
                 arg7_l);
    atomicAdd(&ind_arg3[0+map4idx*18],arg7_l[0]);
    atomicAdd(&ind_arg3[1+map4idx*18],arg7_l[1]);
    atomicAdd(&ind_arg3[2+map4idx*18],arg7_l[2]);
    atomicAdd(&ind_arg3[3+map4idx*18],arg7_l[3]);
    atomicAdd(&ind_arg3[4+map4idx*18],arg7_l[4]);
    atomicAdd(&ind_arg3[5+map4idx*18],arg7_l[5]);
    atomicAdd(&ind_arg3[6+map4idx*18],arg7_l[6]);
    atomicAdd(&ind_arg3[7+map4idx*18],arg7_l[7]);
    atomicAdd(&ind_arg3[8+map4idx*18],arg7_l[8]);
    atomicAdd(&ind_arg3[9+map4idx*18],arg7_l[9]);
    atomicAdd(&ind_arg3[10+map4idx*18],arg7_l[10]);
    atomicAdd(&ind_arg3[11+map4idx*18],arg7_l[11]);
    atomicAdd(&ind_arg3[12+map4idx*18],arg7_l[12]);
    atomicAdd(&ind_arg3[13+map4idx*18],arg7_l[13]);
    atomicAdd(&ind_arg3[14+map4idx*18],arg7_l[14]);
    atomicAdd(&ind_arg3[15+map4idx*18],arg7_l[15]);
    atomicAdd(&ind_arg3[16+map4idx*18],arg7_l[16]);
    atomicAdd(&ind_arg3[17+map4idx*18],arg7_l[17]);
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
  op_arg arg6,
  op_arg arg7){

  double*arg2h = (double *)arg2.data;
  int*arg3h = (int *)arg3.data;
  int nargs = 8;
  op_arg args[8];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(37);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[37].name      = name;
  OP_kernels[37].count    += 1;


  int    ninds   = 4;
  int    inds[8] = {-1,-1,-1,-1,0,1,2,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_bc2\n");
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
    #ifdef OP_BLOCK_SIZE_37
      int nthread = OP_BLOCK_SIZE_37;
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
        op_cuda_pressure_bc2<<<nblocks,nthread>>>(
        (double *)arg4.data_d,
        (double *)arg5.data_d,
        (double *)arg6.data_d,
        (double *)arg7.data_d,
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
  OP_kernels[37].time     += wall_t2 - wall_t1;
}
