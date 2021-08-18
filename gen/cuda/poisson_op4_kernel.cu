//
// auto-generated by op2.py
//

//user function
__device__ void poisson_op4_gpu( const double *cJ, const double *factor, double *op, double *tmp) {
  double cTmp[36 * 10];
  double mm[10 * 10];
  for(int m = 0; m < 36; m++) {
    for(int n = 0; n < 10; n++) {
      int ind = m * 10 + n;
      cTmp[ind] = factor[m] * cJ[m] * cubW_g_cuda[m] * cubV_g_cuda[ind];
    }
  }

  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      mm[c_ind] = 0.0;
      for(int k = 0; k < 36; k++) {
        int b_ind = k * 10 + j;

        int ind = i * 36 + k;
        int a_ind = ((ind * 10) % (10 * 36)) + (ind / 36);

        mm[c_ind] += cubV_g_cuda[b_ind] * cTmp[a_ind];
      }
    }
  }

  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      op[c_ind] += mm[c_ind];
      tmp[c_ind] = mm[c_ind];



    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_op4(
  const double *__restrict arg0,
  const double *__restrict arg1,
  double *arg2,
  double *arg3,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    poisson_op4_gpu(arg0+n*36,
                arg1+n*36,
                arg2+n*100,
                arg3+n*100);
  }
}


//host stub function
void op_par_loop_poisson_op4(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(25);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[25].name      = name;
  OP_kernels[25].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_op4");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_25
      int nthread = OP_BLOCK_SIZE_25;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_poisson_op4<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[25].time     += wall_t2 - wall_t1;
  OP_kernels[25].transfer += (float)set->size * arg0.size;
  OP_kernels[25].transfer += (float)set->size * arg1.size;
  OP_kernels[25].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[25].transfer += (float)set->size * arg3.size * 2.0f;
}
