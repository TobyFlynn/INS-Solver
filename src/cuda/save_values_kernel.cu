//
// auto-generated by op2.py
//

//user function
__device__ void save_values_gpu( const double *v_vals, double *c_vals) {
  c_vals[0]  = (v_vals[0] + v_vals[1] + v_vals[5]) / 3.0;
  c_vals[1]  = (v_vals[1] + v_vals[5] + v_vals[6]) / 3.0;
  c_vals[2]  = (v_vals[1] + v_vals[2] + v_vals[6]) / 3.0;
  c_vals[3]  = (v_vals[2] + v_vals[6] + v_vals[7]) / 3.0;
  c_vals[4]  = (v_vals[2] + v_vals[3] + v_vals[7]) / 3.0;
  c_vals[5]  = (v_vals[3] + v_vals[7] + v_vals[8]) / 3.0;
  c_vals[6]  = (v_vals[3] + v_vals[4] + v_vals[8]) / 3.0;
  c_vals[7]  = (v_vals[5] + v_vals[6] + v_vals[9]) / 3.0;
  c_vals[8]  = (v_vals[6] + v_vals[9] + v_vals[10]) / 3.0;
  c_vals[9]  = (v_vals[6] + v_vals[7] + v_vals[10]) / 3.0;
  c_vals[10] = (v_vals[7] + v_vals[10] + v_vals[11]) / 3.0;
  c_vals[11] = (v_vals[7] + v_vals[8] + v_vals[11]) / 3.0;
  c_vals[12] = (v_vals[9] + v_vals[10] + v_vals[12]) / 3.0;
  c_vals[13] = (v_vals[10] + v_vals[12] + v_vals[13]) / 3.0;
  c_vals[14] = (v_vals[10] + v_vals[11] + v_vals[13]) / 3.0;
  c_vals[15] = (v_vals[12] + v_vals[13] + v_vals[14]) / 3.0;

}

// CUDA kernel function
__global__ void op_cuda_save_values(
  const double *__restrict arg0,
  double *arg1,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    save_values_gpu(arg0+n*15,
                arg1+n*16);
  }
}


//host stub function
void op_par_loop_save_values(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(40);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[40].name      = name;
  OP_kernels[40].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  save_values");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_40
      int nthread = OP_BLOCK_SIZE_40;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_save_values<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[40].time     += wall_t2 - wall_t1;
  OP_kernels[40].transfer += (float)set->size * arg0.size;
  OP_kernels[40].transfer += (float)set->size * arg1.size * 2.0f;
}
