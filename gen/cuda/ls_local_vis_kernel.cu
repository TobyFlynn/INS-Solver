//
// auto-generated by op2.py
//

//user function
__device__ void ls_local_vis_gpu( const double *visMax, const double *modal,
                         double *viscosity) {

  double q[DG_ORDER + 1];
  #if DG_ORDER == 4
  q[0] = modal[0];
  q[1] = modal[1] * modal[1] + modal[5] * modal[5];
  q[1] = sqrt(q[1] / 2.0);
  q[2] = modal[2] * modal[2] + modal[6] * modal[6] + modal[9] * modal[9];
  q[2] = sqrt(q[2] / 3.0);
  q[3] = modal[3] * modal[3] + modal[7] * modal[7] + modal[10] * modal[10]
         + modal[12] * modal[12];
  q[3] = sqrt(q[3] / 4.0);
  q[4] = modal[4] * modal[4] + modal[8] * modal[8] + modal[11] * modal[11]
         + modal[13] * modal[13] + modal[14] * modal[14];
  q[4] = sqrt(q[4] / 5.0);
  #elif DG_ORDER == 3
  q[0] = modal[0];
  q[1] = modal[1] * modal[1] + modal[4] * modal[4];
  q[1] = sqrt(q[1] / 2.0);
  q[2] = modal[2] * modal[2] + modal[5] * modal[5] + modal[7] * modal[7];
  q[2] = sqrt(q[2] / 3.0);
  q[3] = modal[3] * modal[3] + modal[6] * modal[6] + modal[8] * modal[8]
         + modal[9] * modal[9];
  q[3] = sqrt(q[3] / 4.0);
  #elif DG_ORDER == 2
  q[0] = modal[0];
  q[1] = modal[1] * modal[1] + modal[3] * modal[3];
  q[1] = sqrt(q[1] / 2.0);
  q[2] = modal[2] * modal[2] + modal[4] * modal[4] + modal[5] * modal[5];
  q[2] = sqrt(q[2] / 3.0);
  #else
  q[0] = modal[0];
  q[1] = modal[1] * modal[1] + modal[2] * modal[2];
  q[1] = sqrt(q[1] / 2.0);
  #endif

  #if DG_ORDER == 4
  q[0] = fabs(q[0]); q[1] = fabs(q[1]); q[2] = fabs(q[2]); q[3] = fabs(q[3]);
  q[4] = fabs(q[4]);
  q[4] = fmax(q[3], q[4]);
  q[3] = fmax(q[3], q[4]);
  q[2] = fmax(q[2], q[3]);
  q[1] = fmax(q[1], q[2]);
  q[0] = fmax(q[0], q[1]);
  #elif DG_ORDER == 3
  q[0] = fabs(q[0]); q[1] = fabs(q[1]); q[2] = fabs(q[2]); q[3] = fabs(q[3]);
  q[3] = fmax(q[2], q[3]);
  q[2] = fmax(q[2], q[3]);
  q[1] = fmax(q[1], q[2]);
  q[0] = fmax(q[0], q[1]);
  #elif DG_ORDER == 2
  q[0] = fabs(q[0]); q[1] = fabs(q[1]); q[2] = fabs(q[2]);
  q[2] = fmax(q[1], q[2]);
  q[1] = fmax(q[1], q[2]);
  q[0] = fmax(q[0], q[1]);
  #else
  q[0] = fabs(q[0]); q[1] = fabs(q[1]);
  q[1] = fmax(q[0], q[1]);
  q[0] = fmax(q[0], q[1]);
  #endif




  #if DG_ORDER == 1
  *viscosity = *visMax;
  return;
  #endif

  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  double sum4 = 0.0;
  for(int i = 1; i < DG_ORDER + 1; i++) {
    double logx = logf(i);
    double logq = logf(q[i]);
    sum1 += logq * logx;
    sum2 += logq;
    sum3 += logx;
    sum4 += logx * logx;
  }
  double b = (DG_ORDER * sum1 - sum2 * sum3) / (DG_ORDER * sum4 - sum3 * sum3);
  double a = (sum2 - b * sum3) / (double)DG_ORDER;
  double decay_exponent = -b;
  const double PI = 3.141592653589793238463;
  if(decay_exponent < 1.0)
    *viscosity = *visMax;
  else if(decay_exponent > 3.0)
    *viscosity = 0.0;
  else
    *viscosity = *visMax * 0.5 * (1.0 + sin(-PI * (decay_exponent - 2.0) / 2.0));

}

// CUDA kernel function
__global__ void op_cuda_ls_local_vis(
  const double *arg0,
  const double *__restrict arg1,
  double *arg2,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    ls_local_vis_gpu(arg0,
                 arg1+n*10,
                 arg2+n*1);
  }
}


//host stub function
void op_par_loop_ls_local_vis(char const *name, op_set set,
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
  op_timing_realloc(38);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[38].name      = name;
  OP_kernels[38].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  ls_local_vis");
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
    #ifdef OP_BLOCK_SIZE_38
      int nthread = OP_BLOCK_SIZE_38;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_ls_local_vis<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[38].time     += wall_t2 - wall_t1;
  OP_kernels[38].transfer += (float)set->size * arg1.size;
  OP_kernels[38].transfer += (float)set->size * arg2.size * 2.0f;
}
