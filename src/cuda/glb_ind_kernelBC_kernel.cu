//
// auto-generated by op2.py
//

//user function
__device__ void glb_ind_kernelBC_gpu( const int *glb, int *glbBC) {
  glbBC[0] = glb[0];

}

// CUDA kernel function
__global__ void op_cuda_glb_ind_kernelBC(
  const int *__restrict ind_arg0,
  const int *__restrict opDat0Map,
  int *arg1,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    int map0idx;
    map0idx = opDat0Map[n + set_size * 0];

    //user-supplied kernel call
    glb_ind_kernelBC_gpu(ind_arg0+map0idx*1,
                     arg1+n*1);
  }
}


//host stub function
void op_par_loop_glb_ind_kernelBC(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(19);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[19].name      = name;
  OP_kernels[19].count    += 1;


  int    ninds   = 1;
  int    inds[2] = {0,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: glb_ind_kernelBC\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_19
      int nthread = OP_BLOCK_SIZE_19;
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
        op_cuda_glb_ind_kernelBC<<<nblocks,nthread>>>(
        (int *)arg0.data_d,
        arg0.map_data_d,
        (int*)arg1.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[19].time     += wall_t2 - wall_t1;
}
