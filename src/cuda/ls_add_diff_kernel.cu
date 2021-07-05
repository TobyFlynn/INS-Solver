//
// auto-generated by op2.py
//

//user function
__device__ void ls_add_diff_gpu( const double *diff, double *rk, double *dsldx,
                        double *dsrdx, double *dsldy, double *dsrdy) {
  for(int i = 0; i < 15; i++) {
    rk[i] = rk[i] + diff[i];
  }

  for(int i = 0; i < 21; i++) {
    dsldx[i] = 0.0;
    dsrdx[i] = 0.0;
    dsldy[i] = 0.0;
    dsrdy[i] = 0.0;
  }

}

// CUDA kernel function
__global__ void op_cuda_ls_add_diff(
  const double *__restrict arg0,
  double *arg1,
  double *arg2,
  double *arg3,
  double *arg4,
  double *arg5,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    ls_add_diff_gpu(arg0+n*15,
                arg1+n*15,
                arg2+n*21,
                arg3+n*21,
                arg4+n*21,
                arg5+n*21);
  }
}


//host stub function
void op_par_loop_ls_add_diff(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(56);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[56].name      = name;
  OP_kernels[56].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_add_diff");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_56
      int nthread = OP_BLOCK_SIZE_56;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_ls_add_diff<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[56].time     += wall_t2 - wall_t1;
  OP_kernels[56].transfer += (float)set->size * arg0.size;
  OP_kernels[56].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[56].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[56].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[56].transfer += (float)set->size * arg4.size * 2.0f;
  OP_kernels[56].transfer += (float)set->size * arg5.size * 2.0f;
}
