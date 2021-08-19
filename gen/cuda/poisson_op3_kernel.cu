//
// auto-generated by op2.py
//

//user function
__device__ void poisson_op3_gpu( const int *edgeType, const int *edgeNum,
                        const int *d0, const int *d1, const int *d2,
                        const double *mD, const double *sJ,
                        const double *h, const double *gFactor,
                        const double *factor, double *op1) {
  if(*edgeType != *d0 && *edgeType != *d1 && *edgeType != *d2)
    return;

  const double *gVM;
  if(*edgeNum == 0) {
    gVM = gFInterp0_g_cuda;
  } else if(*edgeNum == 1) {
    gVM = gFInterp1_g_cuda;
  } else {
    gVM = gFInterp2_g_cuda;
  }


  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      for(int k = 0; k < 6; k++) {

        int b_ind = k * 10 + j;

        int a_ind = k * 10 + i;

        int factors_ind = *edgeNum * 6 + k;



        op1[c_ind] += -gVM[a_ind] * gaussW_g_cuda[k] * sJ[factors_ind]
                      * gFactor[factors_ind] * mD[b_ind];


      }
    }
  }


  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      for(int k = 0; k < 6; k++) {

        int b_ind = k * 10 + j;

        int a_ind = k * 10 + i;

        int factors_ind = *edgeNum * 6 + k;



        op1[c_ind] += -gFactor[factors_ind] * mD[a_ind] * gaussW_g_cuda[k]
                      * sJ[factors_ind] * gVM[b_ind];




      }
    }
  }

  double tauA[6];
  double maxTau = 0.0;
  for(int i = 0; i < 6; i++) {
    int ind = *edgeNum  * 6 + i;

    tauA[i] = (DG_ORDER + 1) * (DG_ORDER + 2) * (*h * gFactor[ind]);


  }






  for(int i = 0; i < 10; i++) {
    for(int j = 0; j < 10; j++) {
      int c_ind = i * 10 + j;
      for(int k = 0; k < 6; k++) {

        int b_ind = k * 10 + j;

        int a_ind = k * 10 + i;

        int factors_ind = *edgeNum * 6 + k;



        op1[c_ind] += gVM[a_ind] * gaussW_g_cuda[k] * sJ[factors_ind]
                      * tauA[k] * gVM[b_ind];
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_op3(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  double *__restrict ind_arg4,
  const int *__restrict opDat6Map,
  const int *__restrict arg0,
  const int *__restrict arg1,
  const int *arg2,
  const int *arg3,
  const int *arg4,
  const double *__restrict arg5,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg10_l[100];
    for ( int d=0; d<100; d++ ){
      arg10_l[d] = ZERO_double;
    }
    int map6idx;
    map6idx = opDat6Map[n + set_size * 0];

    //user-supplied kernel call
    poisson_op3_gpu(arg0+n*1,
                arg1+n*1,
                arg2,
                arg3,
                arg4,
                arg5+n*60,
                ind_arg0+map6idx*18,
                ind_arg1+map6idx*1,
                ind_arg2+map6idx*18,
                ind_arg3+map6idx*10,
                arg10_l);
    atomicAdd(&ind_arg4[0+map6idx*100],arg10_l[0]);
    atomicAdd(&ind_arg4[1+map6idx*100],arg10_l[1]);
    atomicAdd(&ind_arg4[2+map6idx*100],arg10_l[2]);
    atomicAdd(&ind_arg4[3+map6idx*100],arg10_l[3]);
    atomicAdd(&ind_arg4[4+map6idx*100],arg10_l[4]);
    atomicAdd(&ind_arg4[5+map6idx*100],arg10_l[5]);
    atomicAdd(&ind_arg4[6+map6idx*100],arg10_l[6]);
    atomicAdd(&ind_arg4[7+map6idx*100],arg10_l[7]);
    atomicAdd(&ind_arg4[8+map6idx*100],arg10_l[8]);
    atomicAdd(&ind_arg4[9+map6idx*100],arg10_l[9]);
    atomicAdd(&ind_arg4[10+map6idx*100],arg10_l[10]);
    atomicAdd(&ind_arg4[11+map6idx*100],arg10_l[11]);
    atomicAdd(&ind_arg4[12+map6idx*100],arg10_l[12]);
    atomicAdd(&ind_arg4[13+map6idx*100],arg10_l[13]);
    atomicAdd(&ind_arg4[14+map6idx*100],arg10_l[14]);
    atomicAdd(&ind_arg4[15+map6idx*100],arg10_l[15]);
    atomicAdd(&ind_arg4[16+map6idx*100],arg10_l[16]);
    atomicAdd(&ind_arg4[17+map6idx*100],arg10_l[17]);
    atomicAdd(&ind_arg4[18+map6idx*100],arg10_l[18]);
    atomicAdd(&ind_arg4[19+map6idx*100],arg10_l[19]);
    atomicAdd(&ind_arg4[20+map6idx*100],arg10_l[20]);
    atomicAdd(&ind_arg4[21+map6idx*100],arg10_l[21]);
    atomicAdd(&ind_arg4[22+map6idx*100],arg10_l[22]);
    atomicAdd(&ind_arg4[23+map6idx*100],arg10_l[23]);
    atomicAdd(&ind_arg4[24+map6idx*100],arg10_l[24]);
    atomicAdd(&ind_arg4[25+map6idx*100],arg10_l[25]);
    atomicAdd(&ind_arg4[26+map6idx*100],arg10_l[26]);
    atomicAdd(&ind_arg4[27+map6idx*100],arg10_l[27]);
    atomicAdd(&ind_arg4[28+map6idx*100],arg10_l[28]);
    atomicAdd(&ind_arg4[29+map6idx*100],arg10_l[29]);
    atomicAdd(&ind_arg4[30+map6idx*100],arg10_l[30]);
    atomicAdd(&ind_arg4[31+map6idx*100],arg10_l[31]);
    atomicAdd(&ind_arg4[32+map6idx*100],arg10_l[32]);
    atomicAdd(&ind_arg4[33+map6idx*100],arg10_l[33]);
    atomicAdd(&ind_arg4[34+map6idx*100],arg10_l[34]);
    atomicAdd(&ind_arg4[35+map6idx*100],arg10_l[35]);
    atomicAdd(&ind_arg4[36+map6idx*100],arg10_l[36]);
    atomicAdd(&ind_arg4[37+map6idx*100],arg10_l[37]);
    atomicAdd(&ind_arg4[38+map6idx*100],arg10_l[38]);
    atomicAdd(&ind_arg4[39+map6idx*100],arg10_l[39]);
    atomicAdd(&ind_arg4[40+map6idx*100],arg10_l[40]);
    atomicAdd(&ind_arg4[41+map6idx*100],arg10_l[41]);
    atomicAdd(&ind_arg4[42+map6idx*100],arg10_l[42]);
    atomicAdd(&ind_arg4[43+map6idx*100],arg10_l[43]);
    atomicAdd(&ind_arg4[44+map6idx*100],arg10_l[44]);
    atomicAdd(&ind_arg4[45+map6idx*100],arg10_l[45]);
    atomicAdd(&ind_arg4[46+map6idx*100],arg10_l[46]);
    atomicAdd(&ind_arg4[47+map6idx*100],arg10_l[47]);
    atomicAdd(&ind_arg4[48+map6idx*100],arg10_l[48]);
    atomicAdd(&ind_arg4[49+map6idx*100],arg10_l[49]);
    atomicAdd(&ind_arg4[50+map6idx*100],arg10_l[50]);
    atomicAdd(&ind_arg4[51+map6idx*100],arg10_l[51]);
    atomicAdd(&ind_arg4[52+map6idx*100],arg10_l[52]);
    atomicAdd(&ind_arg4[53+map6idx*100],arg10_l[53]);
    atomicAdd(&ind_arg4[54+map6idx*100],arg10_l[54]);
    atomicAdd(&ind_arg4[55+map6idx*100],arg10_l[55]);
    atomicAdd(&ind_arg4[56+map6idx*100],arg10_l[56]);
    atomicAdd(&ind_arg4[57+map6idx*100],arg10_l[57]);
    atomicAdd(&ind_arg4[58+map6idx*100],arg10_l[58]);
    atomicAdd(&ind_arg4[59+map6idx*100],arg10_l[59]);
    atomicAdd(&ind_arg4[60+map6idx*100],arg10_l[60]);
    atomicAdd(&ind_arg4[61+map6idx*100],arg10_l[61]);
    atomicAdd(&ind_arg4[62+map6idx*100],arg10_l[62]);
    atomicAdd(&ind_arg4[63+map6idx*100],arg10_l[63]);
    atomicAdd(&ind_arg4[64+map6idx*100],arg10_l[64]);
    atomicAdd(&ind_arg4[65+map6idx*100],arg10_l[65]);
    atomicAdd(&ind_arg4[66+map6idx*100],arg10_l[66]);
    atomicAdd(&ind_arg4[67+map6idx*100],arg10_l[67]);
    atomicAdd(&ind_arg4[68+map6idx*100],arg10_l[68]);
    atomicAdd(&ind_arg4[69+map6idx*100],arg10_l[69]);
    atomicAdd(&ind_arg4[70+map6idx*100],arg10_l[70]);
    atomicAdd(&ind_arg4[71+map6idx*100],arg10_l[71]);
    atomicAdd(&ind_arg4[72+map6idx*100],arg10_l[72]);
    atomicAdd(&ind_arg4[73+map6idx*100],arg10_l[73]);
    atomicAdd(&ind_arg4[74+map6idx*100],arg10_l[74]);
    atomicAdd(&ind_arg4[75+map6idx*100],arg10_l[75]);
    atomicAdd(&ind_arg4[76+map6idx*100],arg10_l[76]);
    atomicAdd(&ind_arg4[77+map6idx*100],arg10_l[77]);
    atomicAdd(&ind_arg4[78+map6idx*100],arg10_l[78]);
    atomicAdd(&ind_arg4[79+map6idx*100],arg10_l[79]);
    atomicAdd(&ind_arg4[80+map6idx*100],arg10_l[80]);
    atomicAdd(&ind_arg4[81+map6idx*100],arg10_l[81]);
    atomicAdd(&ind_arg4[82+map6idx*100],arg10_l[82]);
    atomicAdd(&ind_arg4[83+map6idx*100],arg10_l[83]);
    atomicAdd(&ind_arg4[84+map6idx*100],arg10_l[84]);
    atomicAdd(&ind_arg4[85+map6idx*100],arg10_l[85]);
    atomicAdd(&ind_arg4[86+map6idx*100],arg10_l[86]);
    atomicAdd(&ind_arg4[87+map6idx*100],arg10_l[87]);
    atomicAdd(&ind_arg4[88+map6idx*100],arg10_l[88]);
    atomicAdd(&ind_arg4[89+map6idx*100],arg10_l[89]);
    atomicAdd(&ind_arg4[90+map6idx*100],arg10_l[90]);
    atomicAdd(&ind_arg4[91+map6idx*100],arg10_l[91]);
    atomicAdd(&ind_arg4[92+map6idx*100],arg10_l[92]);
    atomicAdd(&ind_arg4[93+map6idx*100],arg10_l[93]);
    atomicAdd(&ind_arg4[94+map6idx*100],arg10_l[94]);
    atomicAdd(&ind_arg4[95+map6idx*100],arg10_l[95]);
    atomicAdd(&ind_arg4[96+map6idx*100],arg10_l[96]);
    atomicAdd(&ind_arg4[97+map6idx*100],arg10_l[97]);
    atomicAdd(&ind_arg4[98+map6idx*100],arg10_l[98]);
    atomicAdd(&ind_arg4[99+map6idx*100],arg10_l[99]);
  }
}


//host stub function
void op_par_loop_poisson_op3(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10){

  int*arg2h = (int *)arg2.data;
  int*arg3h = (int *)arg3.data;
  int*arg4h = (int *)arg4.data;
  int nargs = 11;
  op_arg args[11];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;
  args[10] = arg10;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(18);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[18].name      = name;
  OP_kernels[18].count    += 1;


  int    ninds   = 5;
  int    inds[11] = {-1,-1,-1,-1,-1,-1,0,1,2,3,4};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_op3\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(int));
    consts_bytes += ROUND_UP(1*sizeof(int));
    consts_bytes += ROUND_UP(1*sizeof(int));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg2.data   = OP_consts_h + consts_bytes;
    arg2.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg2.data)[d] = arg2h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    arg3.data   = OP_consts_h + consts_bytes;
    arg3.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg3.data)[d] = arg3h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    arg4.data   = OP_consts_h + consts_bytes;
    arg4.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg4.data)[d] = arg4h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_18
      int nthread = OP_BLOCK_SIZE_18;
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
        op_cuda_poisson_op3<<<nblocks,nthread>>>(
        (double *)arg6.data_d,
        (double *)arg7.data_d,
        (double *)arg8.data_d,
        (double *)arg9.data_d,
        (double *)arg10.data_d,
        arg6.map_data_d,
        (int*)arg0.data_d,
        (int*)arg1.data_d,
        (int*)arg2.data_d,
        (int*)arg3.data_d,
        (int*)arg4.data_d,
        (double*)arg5.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[18].time     += wall_t2 - wall_t1;
}
