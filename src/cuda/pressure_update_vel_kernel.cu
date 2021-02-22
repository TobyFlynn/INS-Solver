//
// auto-generated by op2.py
//

//user function
__device__ void pressure_update_vel_gpu( const double *factor, const double *dpdx,
                                const double *dpdy, const double *qt0,
                                const double *qt1, double *qtt0, double *qtt1) {
  for(int i = 0; i < 15; i++) {
    qtt0[i] = qt0[i] - *factor * dpdx[i];
    qtt1[i] = qt1[i] - *factor * dpdy[i];
  }

}

// CUDA kernel function
__global__ void op_cuda_pressure_update_vel(
  const double *arg0,
  const double *__restrict arg1,
  const double *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  double *arg5,
  double *arg6,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    pressure_update_vel_gpu(arg0,
                        arg1+n*15,
                        arg2+n*15,
                        arg3+n*15,
                        arg4+n*15,
                        arg5+n*15,
                        arg6+n*15);
  }
}


//host stub function
void op_par_loop_pressure_update_vel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6){

  double*arg0h = (double *)arg0.data;
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
  op_timing_realloc(12);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[12].name      = name;
  OP_kernels[12].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pressure_update_vel");
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
    #ifdef OP_BLOCK_SIZE_12
      int nthread = OP_BLOCK_SIZE_12;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_pressure_update_vel<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      (double *) arg4.data_d,
      (double *) arg5.data_d,
      (double *) arg6.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[12].time     += wall_t2 - wall_t1;
  OP_kernels[12].transfer += (float)set->size * arg1.size;
  OP_kernels[12].transfer += (float)set->size * arg2.size;
  OP_kernels[12].transfer += (float)set->size * arg3.size;
  OP_kernels[12].transfer += (float)set->size * arg4.size;
  OP_kernels[12].transfer += (float)set->size * arg5.size * 2.0f;
  OP_kernels[12].transfer += (float)set->size * arg6.size * 2.0f;
}
