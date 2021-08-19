//
// auto-generated by op2.py
//

//user function
__device__ void gauss_gfi_faces2_gpu( const int *edgeNum, const bool *rev,
                             double *gVPL, double *gVPR) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;

  const double *gFL, *gFR;
  if(edgeL == 0) {
    gFR = gFInterp0_g_cuda;
  } else if(edgeL == 1) {
    gFR = gFInterp1_g_cuda;
  } else {
    gFR = gFInterp2_g_cuda;
  }

  if(edgeR == 0) {
    gFL = gFInterp0_g_cuda;
  } else if(edgeR == 1) {
    gFL = gFInterp1_g_cuda;
  } else {
    gFL = gFInterp2_g_cuda;
  }

  for(int i = 0; i < 6 * 10; i++) {
    gVPL[i] = 0.0;
    gVPR[i] = 0.0;
  }

  for(int m = 0; m < 6; m++) {
    for(int n = 0; n < 10; n++) {
      int indL, indR;
      if(!reverse) {
        indL = m * 10 + n;
        indR = m * 10 + n;
      } else {
        indL = m * 10 + n;
        indR = (6 - 1 - m) * 10 + n;
      }

      gVPL[indL] += gFL[indR];
      gVPR[indR] += gFR[indL];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_gauss_gfi_faces2(
  const int *__restrict arg0,
  const bool *__restrict arg1,
  double *arg2,
  double *arg3,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    gauss_gfi_faces2_gpu(arg0+n*2,
                     arg1+n*1,
                     arg2+n*60,
                     arg3+n*60);
  }
}


//host stub function
void op_par_loop_gauss_gfi_faces2(char const *name, op_set set,
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
  op_timing_realloc(8);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[8].name      = name;
  OP_kernels[8].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  gauss_gfi_faces2");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_8
      int nthread = OP_BLOCK_SIZE_8;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_gauss_gfi_faces2<<<nblocks,nthread>>>(
      (int *) arg0.data_d,
      (bool *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[8].time     += wall_t2 - wall_t1;
  OP_kernels[8].transfer += (float)set->size * arg0.size;
  OP_kernels[8].transfer += (float)set->size * arg1.size;
  OP_kernels[8].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[8].transfer += (float)set->size * arg3.size * 2.0f;
}
