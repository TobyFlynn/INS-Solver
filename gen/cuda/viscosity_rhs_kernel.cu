//
// auto-generated by op2.py
//

//user function
__device__ void viscosity_rhs_gpu( const double *factor, double *vRHS0, double *vRHS1) {
  for(int i = 0; i < 3; i++) {
    vRHS0[i] = (*factor) * vRHS0[i];
    vRHS1[i] = (*factor) * vRHS1[i];
  }

}

// CUDA kernel function
__global__ void op_cuda_viscosity_rhs(
  const double *arg0,
  double *arg1,
  double *arg2,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    viscosity_rhs_gpu(arg0,
                  arg1+n*3,
                  arg2+n*3);
  }
}


//host stub function
void op_par_loop_viscosity_rhs(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  double*arg0h = (double *)arg0.data;
  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(40);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[40].name      = name;
  OP_kernels[40].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  viscosity_rhs");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
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
    #ifdef OP_BLOCK_SIZE_40
      int nthread = OP_BLOCK_SIZE_40;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_viscosity_rhs<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[40].time     += wall_t2 - wall_t1;
  OP_kernels[40].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[40].transfer += (float)set->size * arg2.size * 2.0f;
}
