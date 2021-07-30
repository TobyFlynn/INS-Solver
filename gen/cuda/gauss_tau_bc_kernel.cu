//
// auto-generated by op2.py
//

//user function
__device__ void gauss_tau_bc_gpu( const int *bedgeNum, const double *fscale, double *tau) {
  tau[*bedgeNum] += 20 * 25 * fscale[*bedgeNum * 2];

}

// CUDA kernel function
__global__ void op_cuda_gauss_tau_bc(
  const double *__restrict ind_arg0,
  double *__restrict ind_arg1,
  const int *__restrict opDat1Map,
  const int *__restrict arg0,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg2_l[3];
    for ( int d=0; d<3; d++ ){
      arg2_l[d] = ZERO_double;
    }
    int map1idx;
    map1idx = opDat1Map[n + set_size * 0];

    //user-supplied kernel call
    gauss_tau_bc_gpu(arg0+n*1,
                 ind_arg0+map1idx*6,
                 arg2_l);
    atomicAdd(&ind_arg1[0+map1idx*3],arg2_l[0]);
    atomicAdd(&ind_arg1[1+map1idx*3],arg2_l[1]);
    atomicAdd(&ind_arg1[2+map1idx*3],arg2_l[2]);
  }
}


//host stub function
void op_par_loop_gauss_tau_bc(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(5);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[5].name      = name;
  OP_kernels[5].count    += 1;


  int    ninds   = 2;
  int    inds[3] = {-1,0,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: gauss_tau_bc\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_5
      int nthread = OP_BLOCK_SIZE_5;
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
        op_cuda_gauss_tau_bc<<<nblocks,nthread>>>(
        (double *)arg1.data_d,
        (double *)arg2.data_d,
        arg1.map_data_d,
        (int*)arg0.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[5].time     += wall_t2 - wall_t1;
}
