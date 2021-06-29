//
// auto-generated by op2.py
//

//user function
__device__ void pressure_solve_1_gpu( const int *edgeNum, const bool *rev,
                             const double *mD0L, const double *mD0R,
                             const double *mD1L, const double *mD1R,
                             const double *mD2L, const double *mD2R,
                             const double *pD0L, const double *pD0R,
                             const double *pD1L, const double *pD1R,
                             const double *pD2L, const double *pD2R,
                             const double *gVP0L, const double *gVP0R,
                             const double *gVP1L, const double *gVP1R,
                             const double *gVP2L, const double *gVP2R,
                             const double *sJL, const double *sJR,
                             const double *hL, const double *hR,
                             const double *tauOL, const double *tauOR,
                             const double *gRhoL, const double *gRhoR,
                             const double *rhoL, const double *rhoR,
                             const double *uL, const double *uR,
                             double *rhsL, double *rhsR) {

  int edgeL = edgeNum[0];
  int edgeR = edgeNum[1];
  bool reverse = *rev;


  const double *mDL, *mDR, *pDL, *pDR, *gVML, *gVMR, *gVPL, *gVPR;
  if(edgeL == 0) {
    mDL  = mD0L;
    pDL  = pD0L;
    gVML = gFInterp0_g_cuda;
    gVPL = gVP0L;
  } else if(edgeL == 1) {
    mDL  = mD1L;
    pDL  = pD1L;
    gVML = gFInterp1_g_cuda;
    gVPL = gVP1L;
  } else {
    mDL  = mD2L;
    pDL  = pD2L;
    gVML = gFInterp2_g_cuda;
    gVPL = gVP2L;
  }

  if(edgeR == 0) {
    mDR  = mD0R;
    pDR  = pD0R;
    gVMR = gFInterp0_g_cuda;
    gVPR = gVP0R;
  } else if(edgeR == 1) {
    mDR  = mD1R;
    pDR  = pD1R;
    gVMR = gFInterp1_g_cuda;
    gVPR = gVP1R;
  } else {
    mDR  = mD2R;
    pDR  = pD2R;
    gVMR = gFInterp2_g_cuda;
    gVPR = gVP2R;
  }

  double op1L[15 * 15];
  double op1R[15 * 15];
  double op2L[15 * 15];
  double op2R[15 * 15];



  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      op1L[c_ind] = 0.0;
      op1R[c_ind] = 0.0;
      op2L[c_ind] = 0.0;
      op2R[c_ind] = 0.0;
      for(int k = 0; k < 7; k++) {

        int b_ind = k * 15 + j;

        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);


        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;
        int factors_indLR;
        int factors_indRR;
        if(reverse) {
          factors_indLR = edgeL * 7 + 6 - k;
          factors_indRR = edgeR * 7 + 6 - k;
        } else {
          factors_indLR = edgeL * 7 + k;
          factors_indRR = edgeR * 7 + k;
        }


        op1L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * (1.0 / gRhoL[factors_indL]) * mDL[b_ind];
        op1R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
                       * (1.0 / gRhoR[factors_indR]) * mDR[b_ind];

        op2L[c_ind] += -0.5 * gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * (1.0 / gRhoR[factors_indRR]) * pDL[b_ind];
        op2R[c_ind] += -0.5 * gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
                       * (1.0 / gRhoL[factors_indLR]) * pDR[b_ind];

      }
    }
  }



  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      for(int k = 0; k < 7; k++) {

        int b_ind = k * 15 + j;

        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);

        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;

        op1L[c_ind] += -(1.0 / rhoL[i]) * mDL[a_ind] * gaussW_g_cuda[k]
                       * sJL[factors_indL] * gVML[b_ind];
        op1R[c_ind] += -(1.0 / rhoR[i]) * mDR[a_ind] * gaussW_g_cuda[k]
                       * sJR[factors_indR] * gVMR[b_ind];

        op2L[c_ind] += (1.0 / rhoL[i]) * mDL[a_ind] * gaussW_g_cuda[k]
                       * sJL[factors_indL] * gVPL[b_ind];
        op2R[c_ind] += (1.0 / rhoR[i]) * mDR[a_ind] * gaussW_g_cuda[k]
                       * sJR[factors_indR] * gVPR[b_ind];

      }
    }
  }

  double tauL[7];
  double tauR[7];
  double maxL = 0.0;
  double maxR = 0.0;
  for(int i = 0; i < 7; i++) {
    int indL = edgeL * 7 + i;
    int indR;
    if(reverse)
      indR = edgeR * 7 + 6 - i;
    else
      indR = edgeR * 7 + i;
    tauL[i] = 10 * 0.5 * 5 * 6 * fmax(*hL / gRhoL[indL], *hR / gRhoR[indR]);
    if(maxL < tauL[i]) {
      maxL = tauL[i];
    }
  }
  for(int i = 0; i < 7; i++) {
    int indL;
    int indR = edgeR * 7 + i;
    if(reverse)
      indL = edgeL * 7 + 6 - i;
    else
      indL = edgeL * 7 + i;
    tauR[i] = 10 * 0.5 * 5 * 6 * fmax(*hL / gRhoL[indL], *hR / gRhoR[indR]);
    if(maxR < tauR[i]) {
      maxR = tauR[i];
    }
  }

  for(int i = 0; i < 7; i++) {
    tauL[i] = maxL;
    tauR[i] = maxR;
  }



  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int c_ind = i * 15 + j;
      for(int k = 0; k < 7; k++) {

        int b_ind = k * 15 + j;

        int ind = i * 7 + k;
        int a_ind = ((ind * 15) % (15 * 7)) + (ind / 7);

        int factors_indL = edgeL * 7 + k;
        int factors_indR = edgeR * 7 + k;

        op1L[c_ind] += gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * tauL[k] * gVML[b_ind];
        op1R[c_ind] += gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
                       * tauR[k] * gVMR[b_ind];

        op2L[c_ind] += -gVML[a_ind] * gaussW_g_cuda[k] * sJL[factors_indL]
                       * tauL[k] * gVPL[b_ind];
        op2R[c_ind] += -gVMR[a_ind] * gaussW_g_cuda[k] * sJR[factors_indR]
                       * tauR[k] * gVPR[b_ind];

      }
    }
  }


  for(int i = 0; i < 15; i++) {
    for(int j = 0; j < 15; j++) {
      int op_ind = i * 15 + j;
      rhsL[i] += op1L[op_ind] * uL[j] + op2L[op_ind] * uR[j];
      rhsR[i] += op1R[op_ind] * uR[j] + op2R[op_ind] * uL[j];
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_pressure_solve_1(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  const double *__restrict ind_arg6,
  const double *__restrict ind_arg7,
  const double *__restrict ind_arg8,
  const double *__restrict ind_arg9,
  const double *__restrict ind_arg10,
  const double *__restrict ind_arg11,
  const double *__restrict ind_arg12,
  const double *__restrict ind_arg13,
  const double *__restrict ind_arg14,
  double *__restrict ind_arg15,
  const int *__restrict opDat2Map,
  const int *__restrict arg0,
  const bool *__restrict arg1,
  int start,
  int end,
  int   set_size) {
  double arg32_l[15];
  double arg33_l[15];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg32_l[15];
    for ( int d=0; d<15; d++ ){
      arg32_l[d] = ZERO_double;
    }
    double arg33_l[15];
    for ( int d=0; d<15; d++ ){
      arg33_l[d] = ZERO_double;
    }
    int map2idx;
    int map3idx;
    map2idx = opDat2Map[n + set_size * 0];
    map3idx = opDat2Map[n + set_size * 1];

    //user-supplied kernel call
    pressure_solve_1_gpu(arg0+n*2,
                     arg1+n*1,
                     ind_arg0+map2idx*105,
                     ind_arg0+map3idx*105,
                     ind_arg1+map2idx*105,
                     ind_arg1+map3idx*105,
                     ind_arg2+map2idx*105,
                     ind_arg2+map3idx*105,
                     ind_arg3+map2idx*105,
                     ind_arg3+map3idx*105,
                     ind_arg4+map2idx*105,
                     ind_arg4+map3idx*105,
                     ind_arg5+map2idx*105,
                     ind_arg5+map3idx*105,
                     ind_arg6+map2idx*105,
                     ind_arg6+map3idx*105,
                     ind_arg7+map2idx*105,
                     ind_arg7+map3idx*105,
                     ind_arg8+map2idx*105,
                     ind_arg8+map3idx*105,
                     ind_arg9+map2idx*21,
                     ind_arg9+map3idx*21,
                     ind_arg10+map2idx*1,
                     ind_arg10+map3idx*1,
                     ind_arg11+map2idx*3,
                     ind_arg11+map3idx*3,
                     ind_arg12+map2idx*21,
                     ind_arg12+map3idx*21,
                     ind_arg13+map2idx*15,
                     ind_arg13+map3idx*15,
                     ind_arg14+map2idx*15,
                     ind_arg14+map3idx*15,
                     arg32_l,
                     arg33_l);
    atomicAdd(&ind_arg15[0+map2idx*15],arg32_l[0]);
    atomicAdd(&ind_arg15[1+map2idx*15],arg32_l[1]);
    atomicAdd(&ind_arg15[2+map2idx*15],arg32_l[2]);
    atomicAdd(&ind_arg15[3+map2idx*15],arg32_l[3]);
    atomicAdd(&ind_arg15[4+map2idx*15],arg32_l[4]);
    atomicAdd(&ind_arg15[5+map2idx*15],arg32_l[5]);
    atomicAdd(&ind_arg15[6+map2idx*15],arg32_l[6]);
    atomicAdd(&ind_arg15[7+map2idx*15],arg32_l[7]);
    atomicAdd(&ind_arg15[8+map2idx*15],arg32_l[8]);
    atomicAdd(&ind_arg15[9+map2idx*15],arg32_l[9]);
    atomicAdd(&ind_arg15[10+map2idx*15],arg32_l[10]);
    atomicAdd(&ind_arg15[11+map2idx*15],arg32_l[11]);
    atomicAdd(&ind_arg15[12+map2idx*15],arg32_l[12]);
    atomicAdd(&ind_arg15[13+map2idx*15],arg32_l[13]);
    atomicAdd(&ind_arg15[14+map2idx*15],arg32_l[14]);
    atomicAdd(&ind_arg15[0+map3idx*15],arg33_l[0]);
    atomicAdd(&ind_arg15[1+map3idx*15],arg33_l[1]);
    atomicAdd(&ind_arg15[2+map3idx*15],arg33_l[2]);
    atomicAdd(&ind_arg15[3+map3idx*15],arg33_l[3]);
    atomicAdd(&ind_arg15[4+map3idx*15],arg33_l[4]);
    atomicAdd(&ind_arg15[5+map3idx*15],arg33_l[5]);
    atomicAdd(&ind_arg15[6+map3idx*15],arg33_l[6]);
    atomicAdd(&ind_arg15[7+map3idx*15],arg33_l[7]);
    atomicAdd(&ind_arg15[8+map3idx*15],arg33_l[8]);
    atomicAdd(&ind_arg15[9+map3idx*15],arg33_l[9]);
    atomicAdd(&ind_arg15[10+map3idx*15],arg33_l[10]);
    atomicAdd(&ind_arg15[11+map3idx*15],arg33_l[11]);
    atomicAdd(&ind_arg15[12+map3idx*15],arg33_l[12]);
    atomicAdd(&ind_arg15[13+map3idx*15],arg33_l[13]);
    atomicAdd(&ind_arg15[14+map3idx*15],arg33_l[14]);
  }
}


//host stub function
void op_par_loop_pressure_solve_1(char const *name, op_set set,
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
  op_arg arg10,
  op_arg arg11,
  op_arg arg12,
  op_arg arg13,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16,
  op_arg arg17,
  op_arg arg18,
  op_arg arg19,
  op_arg arg20,
  op_arg arg21,
  op_arg arg22,
  op_arg arg23,
  op_arg arg24,
  op_arg arg25,
  op_arg arg26,
  op_arg arg27,
  op_arg arg28,
  op_arg arg29,
  op_arg arg30,
  op_arg arg31,
  op_arg arg32,
  op_arg arg33){

  int nargs = 34;
  op_arg args[34];

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
  args[11] = arg11;
  args[12] = arg12;
  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;
  args[16] = arg16;
  args[17] = arg17;
  args[18] = arg18;
  args[19] = arg19;
  args[20] = arg20;
  args[21] = arg21;
  args[22] = arg22;
  args[23] = arg23;
  args[24] = arg24;
  args[25] = arg25;
  args[26] = arg26;
  args[27] = arg27;
  args[28] = arg28;
  args[29] = arg29;
  args[30] = arg30;
  args[31] = arg31;
  args[32] = arg32;
  args[33] = arg33;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(28);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[28].name      = name;
  OP_kernels[28].count    += 1;


  int    ninds   = 16;
  int    inds[34] = {-1,-1,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pressure_solve_1\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_28
      int nthread = OP_BLOCK_SIZE_28;
    #else
      int nthread = OP_block_size;
    #endif

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        int nblocks = (end-start-1)/nthread+1;
        op_cuda_pressure_solve_1<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg4.data_d,
        (double *)arg6.data_d,
        (double *)arg8.data_d,
        (double *)arg10.data_d,
        (double *)arg12.data_d,
        (double *)arg14.data_d,
        (double *)arg16.data_d,
        (double *)arg18.data_d,
        (double *)arg20.data_d,
        (double *)arg22.data_d,
        (double *)arg24.data_d,
        (double *)arg26.data_d,
        (double *)arg28.data_d,
        (double *)arg30.data_d,
        (double *)arg32.data_d,
        arg2.map_data_d,
        (int*)arg0.data_d,
        (bool*)arg1.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[28].time     += wall_t2 - wall_t1;
}
