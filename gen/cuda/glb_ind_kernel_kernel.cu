//
// auto-generated by op2.py
//

//user function
__device__ void glb_ind_kernel_gpu( const int **glb, int *glbL, int *glbR) {
  glbL[0] = glb[0][0];
  glbR[0] = glb[1][0];

}

// CUDA kernel function
__global__ void op_cuda_glb_ind_kernel(
  const int *__restrict ind_arg0,
  const int *__restrict opDat0Map,
  int *arg2,
  int *arg3,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    int map0idx;
    int map1idx;
    map0idx = opDat0Map[n + set_size * 0];
    map1idx = opDat0Map[n + set_size * 1];
    const int* arg0_vec[] = {
       &ind_arg0[1 * map0idx],
       &ind_arg0[1 * map1idx]};

    //user-supplied kernel call
    glb_ind_kernel_gpu(arg0_vec,
                   arg2+n*1,
                   arg3+n*1);
  }
}


//host stub function
void op_par_loop_glb_ind_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  arg0.idx = 0;
  args[0] = arg0;
  for ( int v=1; v<2; v++ ){
    args[0 + v] = op_arg_dat(arg0.dat, v, arg0.map, 1, "int", OP_READ);
  }

  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(15);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[15].name      = name;
  OP_kernels[15].count    += 1;


  int    ninds   = 1;
  int    inds[4] = {0,0,-1,-1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: glb_ind_kernel\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_15
      int nthread = OP_BLOCK_SIZE_15;
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
        op_cuda_glb_ind_kernel<<<nblocks,nthread>>>(
        (int *)arg0.data_d,
        arg0.map_data_d,
        (int*)arg2.data_d,
        (int*)arg3.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[15].time     += wall_t2 - wall_t1;
}