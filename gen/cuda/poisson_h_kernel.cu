//
// auto-generated by op2.py
//

//user function
__device__ void poisson_h_gpu( const double *x, const double *y, double *h) {
  double len0 = sqrt((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1]));
  double len1 = sqrt((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2]));
  double len2 = sqrt((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0]));
  double sper = (len0 + len1 + len2) / 2.0;
  double area = sqrt(sper * (sper - len0) * (sper - len1) * (sper - len2));
  *h = 2.0 * sper / area;

}

// CUDA kernel function
__global__ void op_cuda_poisson_h(
  const double *__restrict arg0,
  const double *__restrict arg1,
  double *arg2,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    poisson_h_gpu(arg0+n*3,
              arg1+n*3,
              arg2+n*1);
  }
}


//host stub function
void op_par_loop_poisson_h(char const *name, op_set set,
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
  op_timing_realloc(41);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[41].name      = name;
  OP_kernels[41].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  poisson_h");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_41
      int nthread = OP_BLOCK_SIZE_41;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_poisson_h<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[41].time     += wall_t2 - wall_t1;
  OP_kernels[41].transfer += (float)set->size * arg0.size;
  OP_kernels[41].transfer += (float)set->size * arg1.size;
  OP_kernels[41].transfer += (float)set->size * arg2.size * 2.0f;
}
