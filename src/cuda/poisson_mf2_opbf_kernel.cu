//
// auto-generated by op2.py
//

//user function
__device__ void poisson_mf2_opbf_gpu( const double *tol, const int *btype, const int *edgeNum,
                             const int *d0, const int *d1, const int *d2, const double *gop0,
                             const double *gop1, const double *gop2, double *op1) {
  if(*btype == *d0 || *btype == *d1 || *btype == *d2) {
    if(*edgeNum == 0) {
      for(int m = 0; m < 15; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          int colInd = n * 15 + m;
          double val = gop0[colInd];
          if(fabs(val) > *tol) {
            op1[ind] += val;
          }
        }
      }
    } else if(*edgeNum == 1) {
      for(int m = 0; m < 15; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          int colInd = n * 15 + m;
          double val = gop1[colInd];
          if(fabs(val) > *tol) {
            op1[ind] += val;
          }
        }
      }
    } else {
      for(int m = 0; m < 15; m++) {
        for(int n = 0; n < 15; n++) {
          int ind = m * 15 + n;
          int colInd = n * 15 + m;
          double val = gop2[colInd];
          if(fabs(val) > *tol) {
            op1[ind] += val;
          }
        }
      }
    }
  }

}

// CUDA kernel function
__global__ void op_cuda_poisson_mf2_opbf(
  const double *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  double *__restrict ind_arg3,
  const int *__restrict opDat6Map,
  const double *arg0,
  const int *__restrict arg1,
  const int *__restrict arg2,
  const int *arg3,
  const int *arg4,
  const int *arg5,
  int start,
  int end,
  int   set_size) {
  double arg9_l[225];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    double arg9_l[225];
    for ( int d=0; d<225; d++ ){
      arg9_l[d] = ZERO_double;
    }
    int map6idx;
    map6idx = opDat6Map[n + set_size * 0];

    //user-supplied kernel call
    poisson_mf2_opbf_gpu(arg0,
                     arg1+n*1,
                     arg2+n*1,
                     arg3,
                     arg4,
                     arg5,
                     ind_arg0+map6idx*225,
                     ind_arg1+map6idx*225,
                     ind_arg2+map6idx*225,
                     arg9_l);
    atomicAdd(&ind_arg3[0+map6idx*225],arg9_l[0]);
    atomicAdd(&ind_arg3[1+map6idx*225],arg9_l[1]);
    atomicAdd(&ind_arg3[2+map6idx*225],arg9_l[2]);
    atomicAdd(&ind_arg3[3+map6idx*225],arg9_l[3]);
    atomicAdd(&ind_arg3[4+map6idx*225],arg9_l[4]);
    atomicAdd(&ind_arg3[5+map6idx*225],arg9_l[5]);
    atomicAdd(&ind_arg3[6+map6idx*225],arg9_l[6]);
    atomicAdd(&ind_arg3[7+map6idx*225],arg9_l[7]);
    atomicAdd(&ind_arg3[8+map6idx*225],arg9_l[8]);
    atomicAdd(&ind_arg3[9+map6idx*225],arg9_l[9]);
    atomicAdd(&ind_arg3[10+map6idx*225],arg9_l[10]);
    atomicAdd(&ind_arg3[11+map6idx*225],arg9_l[11]);
    atomicAdd(&ind_arg3[12+map6idx*225],arg9_l[12]);
    atomicAdd(&ind_arg3[13+map6idx*225],arg9_l[13]);
    atomicAdd(&ind_arg3[14+map6idx*225],arg9_l[14]);
    atomicAdd(&ind_arg3[15+map6idx*225],arg9_l[15]);
    atomicAdd(&ind_arg3[16+map6idx*225],arg9_l[16]);
    atomicAdd(&ind_arg3[17+map6idx*225],arg9_l[17]);
    atomicAdd(&ind_arg3[18+map6idx*225],arg9_l[18]);
    atomicAdd(&ind_arg3[19+map6idx*225],arg9_l[19]);
    atomicAdd(&ind_arg3[20+map6idx*225],arg9_l[20]);
    atomicAdd(&ind_arg3[21+map6idx*225],arg9_l[21]);
    atomicAdd(&ind_arg3[22+map6idx*225],arg9_l[22]);
    atomicAdd(&ind_arg3[23+map6idx*225],arg9_l[23]);
    atomicAdd(&ind_arg3[24+map6idx*225],arg9_l[24]);
    atomicAdd(&ind_arg3[25+map6idx*225],arg9_l[25]);
    atomicAdd(&ind_arg3[26+map6idx*225],arg9_l[26]);
    atomicAdd(&ind_arg3[27+map6idx*225],arg9_l[27]);
    atomicAdd(&ind_arg3[28+map6idx*225],arg9_l[28]);
    atomicAdd(&ind_arg3[29+map6idx*225],arg9_l[29]);
    atomicAdd(&ind_arg3[30+map6idx*225],arg9_l[30]);
    atomicAdd(&ind_arg3[31+map6idx*225],arg9_l[31]);
    atomicAdd(&ind_arg3[32+map6idx*225],arg9_l[32]);
    atomicAdd(&ind_arg3[33+map6idx*225],arg9_l[33]);
    atomicAdd(&ind_arg3[34+map6idx*225],arg9_l[34]);
    atomicAdd(&ind_arg3[35+map6idx*225],arg9_l[35]);
    atomicAdd(&ind_arg3[36+map6idx*225],arg9_l[36]);
    atomicAdd(&ind_arg3[37+map6idx*225],arg9_l[37]);
    atomicAdd(&ind_arg3[38+map6idx*225],arg9_l[38]);
    atomicAdd(&ind_arg3[39+map6idx*225],arg9_l[39]);
    atomicAdd(&ind_arg3[40+map6idx*225],arg9_l[40]);
    atomicAdd(&ind_arg3[41+map6idx*225],arg9_l[41]);
    atomicAdd(&ind_arg3[42+map6idx*225],arg9_l[42]);
    atomicAdd(&ind_arg3[43+map6idx*225],arg9_l[43]);
    atomicAdd(&ind_arg3[44+map6idx*225],arg9_l[44]);
    atomicAdd(&ind_arg3[45+map6idx*225],arg9_l[45]);
    atomicAdd(&ind_arg3[46+map6idx*225],arg9_l[46]);
    atomicAdd(&ind_arg3[47+map6idx*225],arg9_l[47]);
    atomicAdd(&ind_arg3[48+map6idx*225],arg9_l[48]);
    atomicAdd(&ind_arg3[49+map6idx*225],arg9_l[49]);
    atomicAdd(&ind_arg3[50+map6idx*225],arg9_l[50]);
    atomicAdd(&ind_arg3[51+map6idx*225],arg9_l[51]);
    atomicAdd(&ind_arg3[52+map6idx*225],arg9_l[52]);
    atomicAdd(&ind_arg3[53+map6idx*225],arg9_l[53]);
    atomicAdd(&ind_arg3[54+map6idx*225],arg9_l[54]);
    atomicAdd(&ind_arg3[55+map6idx*225],arg9_l[55]);
    atomicAdd(&ind_arg3[56+map6idx*225],arg9_l[56]);
    atomicAdd(&ind_arg3[57+map6idx*225],arg9_l[57]);
    atomicAdd(&ind_arg3[58+map6idx*225],arg9_l[58]);
    atomicAdd(&ind_arg3[59+map6idx*225],arg9_l[59]);
    atomicAdd(&ind_arg3[60+map6idx*225],arg9_l[60]);
    atomicAdd(&ind_arg3[61+map6idx*225],arg9_l[61]);
    atomicAdd(&ind_arg3[62+map6idx*225],arg9_l[62]);
    atomicAdd(&ind_arg3[63+map6idx*225],arg9_l[63]);
    atomicAdd(&ind_arg3[64+map6idx*225],arg9_l[64]);
    atomicAdd(&ind_arg3[65+map6idx*225],arg9_l[65]);
    atomicAdd(&ind_arg3[66+map6idx*225],arg9_l[66]);
    atomicAdd(&ind_arg3[67+map6idx*225],arg9_l[67]);
    atomicAdd(&ind_arg3[68+map6idx*225],arg9_l[68]);
    atomicAdd(&ind_arg3[69+map6idx*225],arg9_l[69]);
    atomicAdd(&ind_arg3[70+map6idx*225],arg9_l[70]);
    atomicAdd(&ind_arg3[71+map6idx*225],arg9_l[71]);
    atomicAdd(&ind_arg3[72+map6idx*225],arg9_l[72]);
    atomicAdd(&ind_arg3[73+map6idx*225],arg9_l[73]);
    atomicAdd(&ind_arg3[74+map6idx*225],arg9_l[74]);
    atomicAdd(&ind_arg3[75+map6idx*225],arg9_l[75]);
    atomicAdd(&ind_arg3[76+map6idx*225],arg9_l[76]);
    atomicAdd(&ind_arg3[77+map6idx*225],arg9_l[77]);
    atomicAdd(&ind_arg3[78+map6idx*225],arg9_l[78]);
    atomicAdd(&ind_arg3[79+map6idx*225],arg9_l[79]);
    atomicAdd(&ind_arg3[80+map6idx*225],arg9_l[80]);
    atomicAdd(&ind_arg3[81+map6idx*225],arg9_l[81]);
    atomicAdd(&ind_arg3[82+map6idx*225],arg9_l[82]);
    atomicAdd(&ind_arg3[83+map6idx*225],arg9_l[83]);
    atomicAdd(&ind_arg3[84+map6idx*225],arg9_l[84]);
    atomicAdd(&ind_arg3[85+map6idx*225],arg9_l[85]);
    atomicAdd(&ind_arg3[86+map6idx*225],arg9_l[86]);
    atomicAdd(&ind_arg3[87+map6idx*225],arg9_l[87]);
    atomicAdd(&ind_arg3[88+map6idx*225],arg9_l[88]);
    atomicAdd(&ind_arg3[89+map6idx*225],arg9_l[89]);
    atomicAdd(&ind_arg3[90+map6idx*225],arg9_l[90]);
    atomicAdd(&ind_arg3[91+map6idx*225],arg9_l[91]);
    atomicAdd(&ind_arg3[92+map6idx*225],arg9_l[92]);
    atomicAdd(&ind_arg3[93+map6idx*225],arg9_l[93]);
    atomicAdd(&ind_arg3[94+map6idx*225],arg9_l[94]);
    atomicAdd(&ind_arg3[95+map6idx*225],arg9_l[95]);
    atomicAdd(&ind_arg3[96+map6idx*225],arg9_l[96]);
    atomicAdd(&ind_arg3[97+map6idx*225],arg9_l[97]);
    atomicAdd(&ind_arg3[98+map6idx*225],arg9_l[98]);
    atomicAdd(&ind_arg3[99+map6idx*225],arg9_l[99]);
    atomicAdd(&ind_arg3[100+map6idx*225],arg9_l[100]);
    atomicAdd(&ind_arg3[101+map6idx*225],arg9_l[101]);
    atomicAdd(&ind_arg3[102+map6idx*225],arg9_l[102]);
    atomicAdd(&ind_arg3[103+map6idx*225],arg9_l[103]);
    atomicAdd(&ind_arg3[104+map6idx*225],arg9_l[104]);
    atomicAdd(&ind_arg3[105+map6idx*225],arg9_l[105]);
    atomicAdd(&ind_arg3[106+map6idx*225],arg9_l[106]);
    atomicAdd(&ind_arg3[107+map6idx*225],arg9_l[107]);
    atomicAdd(&ind_arg3[108+map6idx*225],arg9_l[108]);
    atomicAdd(&ind_arg3[109+map6idx*225],arg9_l[109]);
    atomicAdd(&ind_arg3[110+map6idx*225],arg9_l[110]);
    atomicAdd(&ind_arg3[111+map6idx*225],arg9_l[111]);
    atomicAdd(&ind_arg3[112+map6idx*225],arg9_l[112]);
    atomicAdd(&ind_arg3[113+map6idx*225],arg9_l[113]);
    atomicAdd(&ind_arg3[114+map6idx*225],arg9_l[114]);
    atomicAdd(&ind_arg3[115+map6idx*225],arg9_l[115]);
    atomicAdd(&ind_arg3[116+map6idx*225],arg9_l[116]);
    atomicAdd(&ind_arg3[117+map6idx*225],arg9_l[117]);
    atomicAdd(&ind_arg3[118+map6idx*225],arg9_l[118]);
    atomicAdd(&ind_arg3[119+map6idx*225],arg9_l[119]);
    atomicAdd(&ind_arg3[120+map6idx*225],arg9_l[120]);
    atomicAdd(&ind_arg3[121+map6idx*225],arg9_l[121]);
    atomicAdd(&ind_arg3[122+map6idx*225],arg9_l[122]);
    atomicAdd(&ind_arg3[123+map6idx*225],arg9_l[123]);
    atomicAdd(&ind_arg3[124+map6idx*225],arg9_l[124]);
    atomicAdd(&ind_arg3[125+map6idx*225],arg9_l[125]);
    atomicAdd(&ind_arg3[126+map6idx*225],arg9_l[126]);
    atomicAdd(&ind_arg3[127+map6idx*225],arg9_l[127]);
    atomicAdd(&ind_arg3[128+map6idx*225],arg9_l[128]);
    atomicAdd(&ind_arg3[129+map6idx*225],arg9_l[129]);
    atomicAdd(&ind_arg3[130+map6idx*225],arg9_l[130]);
    atomicAdd(&ind_arg3[131+map6idx*225],arg9_l[131]);
    atomicAdd(&ind_arg3[132+map6idx*225],arg9_l[132]);
    atomicAdd(&ind_arg3[133+map6idx*225],arg9_l[133]);
    atomicAdd(&ind_arg3[134+map6idx*225],arg9_l[134]);
    atomicAdd(&ind_arg3[135+map6idx*225],arg9_l[135]);
    atomicAdd(&ind_arg3[136+map6idx*225],arg9_l[136]);
    atomicAdd(&ind_arg3[137+map6idx*225],arg9_l[137]);
    atomicAdd(&ind_arg3[138+map6idx*225],arg9_l[138]);
    atomicAdd(&ind_arg3[139+map6idx*225],arg9_l[139]);
    atomicAdd(&ind_arg3[140+map6idx*225],arg9_l[140]);
    atomicAdd(&ind_arg3[141+map6idx*225],arg9_l[141]);
    atomicAdd(&ind_arg3[142+map6idx*225],arg9_l[142]);
    atomicAdd(&ind_arg3[143+map6idx*225],arg9_l[143]);
    atomicAdd(&ind_arg3[144+map6idx*225],arg9_l[144]);
    atomicAdd(&ind_arg3[145+map6idx*225],arg9_l[145]);
    atomicAdd(&ind_arg3[146+map6idx*225],arg9_l[146]);
    atomicAdd(&ind_arg3[147+map6idx*225],arg9_l[147]);
    atomicAdd(&ind_arg3[148+map6idx*225],arg9_l[148]);
    atomicAdd(&ind_arg3[149+map6idx*225],arg9_l[149]);
    atomicAdd(&ind_arg3[150+map6idx*225],arg9_l[150]);
    atomicAdd(&ind_arg3[151+map6idx*225],arg9_l[151]);
    atomicAdd(&ind_arg3[152+map6idx*225],arg9_l[152]);
    atomicAdd(&ind_arg3[153+map6idx*225],arg9_l[153]);
    atomicAdd(&ind_arg3[154+map6idx*225],arg9_l[154]);
    atomicAdd(&ind_arg3[155+map6idx*225],arg9_l[155]);
    atomicAdd(&ind_arg3[156+map6idx*225],arg9_l[156]);
    atomicAdd(&ind_arg3[157+map6idx*225],arg9_l[157]);
    atomicAdd(&ind_arg3[158+map6idx*225],arg9_l[158]);
    atomicAdd(&ind_arg3[159+map6idx*225],arg9_l[159]);
    atomicAdd(&ind_arg3[160+map6idx*225],arg9_l[160]);
    atomicAdd(&ind_arg3[161+map6idx*225],arg9_l[161]);
    atomicAdd(&ind_arg3[162+map6idx*225],arg9_l[162]);
    atomicAdd(&ind_arg3[163+map6idx*225],arg9_l[163]);
    atomicAdd(&ind_arg3[164+map6idx*225],arg9_l[164]);
    atomicAdd(&ind_arg3[165+map6idx*225],arg9_l[165]);
    atomicAdd(&ind_arg3[166+map6idx*225],arg9_l[166]);
    atomicAdd(&ind_arg3[167+map6idx*225],arg9_l[167]);
    atomicAdd(&ind_arg3[168+map6idx*225],arg9_l[168]);
    atomicAdd(&ind_arg3[169+map6idx*225],arg9_l[169]);
    atomicAdd(&ind_arg3[170+map6idx*225],arg9_l[170]);
    atomicAdd(&ind_arg3[171+map6idx*225],arg9_l[171]);
    atomicAdd(&ind_arg3[172+map6idx*225],arg9_l[172]);
    atomicAdd(&ind_arg3[173+map6idx*225],arg9_l[173]);
    atomicAdd(&ind_arg3[174+map6idx*225],arg9_l[174]);
    atomicAdd(&ind_arg3[175+map6idx*225],arg9_l[175]);
    atomicAdd(&ind_arg3[176+map6idx*225],arg9_l[176]);
    atomicAdd(&ind_arg3[177+map6idx*225],arg9_l[177]);
    atomicAdd(&ind_arg3[178+map6idx*225],arg9_l[178]);
    atomicAdd(&ind_arg3[179+map6idx*225],arg9_l[179]);
    atomicAdd(&ind_arg3[180+map6idx*225],arg9_l[180]);
    atomicAdd(&ind_arg3[181+map6idx*225],arg9_l[181]);
    atomicAdd(&ind_arg3[182+map6idx*225],arg9_l[182]);
    atomicAdd(&ind_arg3[183+map6idx*225],arg9_l[183]);
    atomicAdd(&ind_arg3[184+map6idx*225],arg9_l[184]);
    atomicAdd(&ind_arg3[185+map6idx*225],arg9_l[185]);
    atomicAdd(&ind_arg3[186+map6idx*225],arg9_l[186]);
    atomicAdd(&ind_arg3[187+map6idx*225],arg9_l[187]);
    atomicAdd(&ind_arg3[188+map6idx*225],arg9_l[188]);
    atomicAdd(&ind_arg3[189+map6idx*225],arg9_l[189]);
    atomicAdd(&ind_arg3[190+map6idx*225],arg9_l[190]);
    atomicAdd(&ind_arg3[191+map6idx*225],arg9_l[191]);
    atomicAdd(&ind_arg3[192+map6idx*225],arg9_l[192]);
    atomicAdd(&ind_arg3[193+map6idx*225],arg9_l[193]);
    atomicAdd(&ind_arg3[194+map6idx*225],arg9_l[194]);
    atomicAdd(&ind_arg3[195+map6idx*225],arg9_l[195]);
    atomicAdd(&ind_arg3[196+map6idx*225],arg9_l[196]);
    atomicAdd(&ind_arg3[197+map6idx*225],arg9_l[197]);
    atomicAdd(&ind_arg3[198+map6idx*225],arg9_l[198]);
    atomicAdd(&ind_arg3[199+map6idx*225],arg9_l[199]);
    atomicAdd(&ind_arg3[200+map6idx*225],arg9_l[200]);
    atomicAdd(&ind_arg3[201+map6idx*225],arg9_l[201]);
    atomicAdd(&ind_arg3[202+map6idx*225],arg9_l[202]);
    atomicAdd(&ind_arg3[203+map6idx*225],arg9_l[203]);
    atomicAdd(&ind_arg3[204+map6idx*225],arg9_l[204]);
    atomicAdd(&ind_arg3[205+map6idx*225],arg9_l[205]);
    atomicAdd(&ind_arg3[206+map6idx*225],arg9_l[206]);
    atomicAdd(&ind_arg3[207+map6idx*225],arg9_l[207]);
    atomicAdd(&ind_arg3[208+map6idx*225],arg9_l[208]);
    atomicAdd(&ind_arg3[209+map6idx*225],arg9_l[209]);
    atomicAdd(&ind_arg3[210+map6idx*225],arg9_l[210]);
    atomicAdd(&ind_arg3[211+map6idx*225],arg9_l[211]);
    atomicAdd(&ind_arg3[212+map6idx*225],arg9_l[212]);
    atomicAdd(&ind_arg3[213+map6idx*225],arg9_l[213]);
    atomicAdd(&ind_arg3[214+map6idx*225],arg9_l[214]);
    atomicAdd(&ind_arg3[215+map6idx*225],arg9_l[215]);
    atomicAdd(&ind_arg3[216+map6idx*225],arg9_l[216]);
    atomicAdd(&ind_arg3[217+map6idx*225],arg9_l[217]);
    atomicAdd(&ind_arg3[218+map6idx*225],arg9_l[218]);
    atomicAdd(&ind_arg3[219+map6idx*225],arg9_l[219]);
    atomicAdd(&ind_arg3[220+map6idx*225],arg9_l[220]);
    atomicAdd(&ind_arg3[221+map6idx*225],arg9_l[221]);
    atomicAdd(&ind_arg3[222+map6idx*225],arg9_l[222]);
    atomicAdd(&ind_arg3[223+map6idx*225],arg9_l[223]);
    atomicAdd(&ind_arg3[224+map6idx*225],arg9_l[224]);
  }
}


//host stub function
void op_par_loop_poisson_mf2_opbf(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  double*arg0h = (double *)arg0.data;
  int*arg3h = (int *)arg3.data;
  int*arg4h = (int *)arg4.data;
  int*arg5h = (int *)arg5.data;
  int nargs = 10;
  op_arg args[10];

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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(39);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[39].name      = name;
  OP_kernels[39].count    += 1;


  int    ninds   = 4;
  int    inds[10] = {-1,-1,-1,-1,-1,-1,0,1,2,3};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: poisson_mf2_opbf\n");
  }
  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    consts_bytes += ROUND_UP(1*sizeof(int));
    consts_bytes += ROUND_UP(1*sizeof(int));
    consts_bytes += ROUND_UP(1*sizeof(int));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg0.data   = OP_consts_h + consts_bytes;
    arg0.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg0.data)[d] = arg0h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
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
    arg5.data   = OP_consts_h + consts_bytes;
    arg5.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((int *)arg5.data)[d] = arg5h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_39
      int nthread = OP_BLOCK_SIZE_39;
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
        op_cuda_poisson_mf2_opbf<<<nblocks,nthread>>>(
        (double *)arg6.data_d,
        (double *)arg7.data_d,
        (double *)arg8.data_d,
        (double *)arg9.data_d,
        arg6.map_data_d,
        (double*)arg0.data_d,
        (int*)arg1.data_d,
        (int*)arg2.data_d,
        (int*)arg3.data_d,
        (int*)arg4.data_d,
        (int*)arg5.data_d,
        start,end,set->size+set->exec_size);
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[39].time     += wall_t2 - wall_t1;
}
