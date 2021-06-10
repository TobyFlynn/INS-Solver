//
// auto-generated by op2.py
//

//user function
__device__ void ls_reinit_check_gpu( const double *alpha, const double *s,
                            const double *dsdx, const double *dsdy,
                            double *res, int *count) {
  for(int i = 0; i < 15; i++) {
    if(fabs(s[i]) < (*alpha)) {
      *res += dsdx[i] * dsdx[i] + dsdy[i] * dsdy[i];
      *count += 1;
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_ls_reinit_check(
  const double *arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  double *arg4,
  int *arg5,
  int   set_size ) {

  double arg4_l[1];
  for ( int d=0; d<1; d++ ){
    arg4_l[d]=ZERO_double;
  }
  int arg5_l[1];
  for ( int d=0; d<1; d++ ){
    arg5_l[d]=ZERO_int;
  }

  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    ls_reinit_check_gpu(arg0,
                    arg1+n*15,
                    arg2+n*15,
                    arg3+n*15,
                    arg4_l,
                    arg5_l);
  }

  //global reductions

  for ( int d=0; d<1; d++ ){
    op_reduction<OP_INC>(&arg4[d+blockIdx.x*1],arg4_l[d]);
  }
  for ( int d=0; d<1; d++ ){
    op_reduction<OP_INC>(&arg5[d+blockIdx.x*1],arg5_l[d]);
  }
}


//host stub function
void op_par_loop_ls_reinit_check(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  double*arg0h = (double *)arg0.data;
  double*arg4h = (double *)arg4.data;
  int*arg5h = (int *)arg5.data;
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
  op_timing_realloc(72);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[72].name      = name;
  OP_kernels[72].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_reinit_check");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg0.data   = OP_consts_h + consts_bytes;
    arg0.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg0.data)[d] = arg0h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_72
      int nthread = OP_BLOCK_SIZE_72;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    //transfer global reduction data to GPU
    int maxblocks = nblocks;
    int reduct_bytes = 0;
    int reduct_size  = 0;
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    reduct_size   = MAX(reduct_size,sizeof(double));
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(int));
    reduct_size   = MAX(reduct_size,sizeof(int));
    reallocReductArrays(reduct_bytes);
    reduct_bytes = 0;
    arg4.data   = OP_reduct_h + reduct_bytes;
    arg4.data_d = OP_reduct_d + reduct_bytes;
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        ((double *)arg4.data)[d+b*1] = ZERO_double;
      }
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    arg5.data   = OP_reduct_h + reduct_bytes;
    arg5.data_d = OP_reduct_d + reduct_bytes;
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        ((int *)arg5.data)[d+b*1] = ZERO_int;
      }
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(int));
    mvReductArraysToDevice(reduct_bytes);

    int nshared = reduct_size*nthread;
    op_cuda_ls_reinit_check<<<nblocks,nthread,nshared>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (int *) arg5.data_d,
      set->size );
    //transfer global reduction data back to CPU
    mvReductArraysToHost(reduct_bytes);
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        arg4h[d] = arg4h[d] + ((double *)arg4.data)[d+b*1];
      }
    }
    arg4.data = (char *)arg4h;
    op_mpi_reduce(&arg4,arg4h);
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        arg5h[d] = arg5h[d] + ((int *)arg5.data)[d+b*1];
      }
    }
    arg5.data = (char *)arg5h;
    op_mpi_reduce(&arg5,arg5h);
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[72].time     += wall_t2 - wall_t1;
  OP_kernels[72].transfer += (float)set->size * arg1.size;
  OP_kernels[72].transfer += (float)set->size * arg2.size;
  OP_kernels[72].transfer += (float)set->size * arg3.size;
}
