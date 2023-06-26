__constant__ int opDat0_custom_pmf_stride_OP2CONSTANT;
int opDat0_custom_pmf_stride_OP2HOST=-1;
__constant__ int opDat10_custom_pmf_stride_OP2CONSTANT;
int opDat10_custom_pmf_stride_OP2HOST=-1;
__constant__ int direct_custom_pmf_stride_OP2CONSTANT;
int direct_custom_pmf_stride_OP2HOST=-1;

template<int p, int dg_np, int dg_npf>
__device__ void _pmf_3d_mult_complete_flux_gpu(const int *fmash_sh, const int *faceNums,
                            const int *fmaskF, const double *nx, const double *ny,
                            const double *nz, const double *fscale,
                            const double *sJ, const double *geof, const double *in,
                            const double **in_p, const double *in_x,
                            const double *in_y, const double *in_z,
                            const double **in_x_p, const double **in_y_p,
                            const double **in_z_p, double *out) {
  #pragma unroll
  for(int i = 0; i < dg_np; i++) {
    out[(i)*opDat0_custom_pmf_stride_OP2CONSTANT] = 0.0;
  }

  double tmp_4[dg_np], tmp_5[dg_np], tmp_6[dg_np];
  const double cell_J = geof[(J_IND)*opDat0_custom_pmf_stride_OP2CONSTANT];
  const double *mass_mat = &dg_Mass_kernel[(p - 1) * 20 * 20];
  for(int i = 0; i < dg_np; i++) {
    double tmp_gemv[3] = {};
    #pragma unroll
    for(int j = 0; j < dg_np; j++) {
      int ind = DG_MAT_IND(i,j, dg_np, dg_np);
      tmp_gemv[0] += mass_mat[ind] * in_x[(j)*opDat0_custom_pmf_stride_OP2CONSTANT];
      tmp_gemv[1] += mass_mat[ind] * in_y[(j)*opDat0_custom_pmf_stride_OP2CONSTANT];
      tmp_gemv[2] += mass_mat[ind] * in_z[(j)*opDat0_custom_pmf_stride_OP2CONSTANT];
    }
     tmp_4[i] = cell_J * tmp_gemv[0];
     tmp_5[i] = cell_J * tmp_gemv[1];
     tmp_6[i] = cell_J * tmp_gemv[2];
  }

  const double *emat_mat = &dg_Emat_kernel[(p - 1) * 4 * 10 * 20];
  for(int i = 0; i < 4; i++) {
    const double gtau = 2.0 * (p + 1) * (p + 2) * fmax(fscale[(i * 2)*direct_custom_pmf_stride_OP2CONSTANT], fscale[(i * 2 + 1)*direct_custom_pmf_stride_OP2CONSTANT]);

    const int findL = faceNums[(2 * i)*direct_custom_pmf_stride_OP2CONSTANT] * dg_npf;
    const int *fmaskL = &fmash_sh[faceNums[(2 * i)*direct_custom_pmf_stride_OP2CONSTANT] * dg_npf];

    for(int j = 0; j < dg_npf; j++) {
      const int fmaskIndL = fmaskL[j];
      const int fmaskIndR = fmaskF[(i * dg_npf + j)*direct_custom_pmf_stride_OP2CONSTANT];
      const double diffL_u = in[(fmaskIndL)*opDat0_custom_pmf_stride_OP2CONSTANT] - in_p[i][(fmaskIndR)*opDat10_custom_pmf_stride_OP2CONSTANT];
      const double diffL_u_x = nx[(i)*direct_custom_pmf_stride_OP2CONSTANT] * (in_x_p[i][(fmaskIndR)*opDat10_custom_pmf_stride_OP2CONSTANT] + in_x[(fmaskIndL)*opDat0_custom_pmf_stride_OP2CONSTANT]);
      const double diffL_u_y = ny[(i)*direct_custom_pmf_stride_OP2CONSTANT] * (in_y_p[i][(fmaskIndR)*opDat10_custom_pmf_stride_OP2CONSTANT] + in_y[(fmaskIndL)*opDat0_custom_pmf_stride_OP2CONSTANT]);
      const double diffL_u_z = nz[(i)*direct_custom_pmf_stride_OP2CONSTANT] * (in_z_p[i][(fmaskIndR)*opDat10_custom_pmf_stride_OP2CONSTANT] + in_z[(fmaskIndL)*opDat0_custom_pmf_stride_OP2CONSTANT]);
      const double diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

      const int indL = findL + j;
      const double int_fact = 0.5 * sJ[(i)*direct_custom_pmf_stride_OP2CONSTANT];
      const double tmp_0 = int_fact * (gtau * diffL_u - diffL_u_grad);
      const double l_tmpL = int_fact * -diffL_u;
      const double tmp_1 = nx[(i)*direct_custom_pmf_stride_OP2CONSTANT] * l_tmpL;
      const double tmp_2 = ny[(i)*direct_custom_pmf_stride_OP2CONSTANT] * l_tmpL;
      const double tmp_3 = nz[(i)*direct_custom_pmf_stride_OP2CONSTANT] * l_tmpL;
      #pragma unroll
      for(int i = 0; i < dg_np; i++) {
        int ind = DG_MAT_IND(i, indL, dg_np, 4 * dg_npf);
        tmp_4[i] += emat_mat[ind] * tmp_1;
        tmp_5[i] += emat_mat[ind] * tmp_2;
        tmp_6[i] += emat_mat[ind] * tmp_3;
        out[(i)*opDat0_custom_pmf_stride_OP2CONSTANT] += emat_mat[ind] * tmp_0;
      }
    }
  }

  #pragma unroll
  for(int n = 0; n < dg_np; n++) {
    const double tmp_x = tmp_4[n];
    const double tmp_y = tmp_5[n];
    const double tmp_z = tmp_6[n];
    tmp_4[n] = geof[(RX_IND)*opDat0_custom_pmf_stride_OP2CONSTANT] * tmp_x + geof[(RY_IND)*opDat0_custom_pmf_stride_OP2CONSTANT] * tmp_y + geof[(RZ_IND)*opDat0_custom_pmf_stride_OP2CONSTANT] * tmp_z;
    tmp_5[n] = geof[(SX_IND)*opDat0_custom_pmf_stride_OP2CONSTANT] * tmp_x + geof[(SY_IND)*opDat0_custom_pmf_stride_OP2CONSTANT] * tmp_y + geof[(SZ_IND)*opDat0_custom_pmf_stride_OP2CONSTANT] * tmp_z;
    tmp_6[n] = geof[(TX_IND)*opDat0_custom_pmf_stride_OP2CONSTANT] * tmp_x + geof[(TY_IND)*opDat0_custom_pmf_stride_OP2CONSTANT] * tmp_y + geof[(TZ_IND)*opDat0_custom_pmf_stride_OP2CONSTANT] * tmp_z;
  }

  const double *dr_mat = &dg_Dr_kernel[(p - 1) * 20 * 20];
  const double *ds_mat = &dg_Ds_kernel[(p - 1) * 20 * 20];
  const double *dt_mat = &dg_Dt_kernel[(p - 1) * 20 * 20];
  for(int i = 0; i < dg_np; i++) {
    double tmp_gemv = 0.0;
    #pragma unroll
    for(int j = 0; j < dg_np; j++) {
      int ind = DG_MAT_IND(j,i, dg_np, dg_np);
      tmp_gemv += dr_mat[ind] * tmp_4[j];
      tmp_gemv += ds_mat[ind] * tmp_5[j];
      tmp_gemv += dt_mat[ind] * tmp_6[j];
    }
     out[(i)*opDat0_custom_pmf_stride_OP2CONSTANT] += tmp_gemv;
  }
}

template<int p>
__global__ void _op_cuda_pmf_3d_mult_complete_flux(
  const int *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  const double *__restrict ind_arg6,
  const double *__restrict ind_arg7,
  const double *__restrict ind_arg8,
  const double *__restrict ind_arg9,
  double *__restrict ind_arg10,
  const int *__restrict opDat0Map,
  const int *__restrict opDat10Map,
  const int *__restrict arg1,
  const int *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  const double *__restrict arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  int start,
  int end,
  int   set_size) {
  const int np = (p + 1) * (p + 2) * (p + 3) / 6;
  const int npf = (p + 1) * (p + 2) / 2;
  __shared__ int fmask_sh[4 * npf];
  for(int i = threadIdx.x; i < 4 * npf; i++) {
    fmask_sh[i] = FMASK_cuda[(p - 1) * 4 * 10 + i];
  }
  __syncthreads();

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    int map0idx;
    int map10idx;
    int map11idx;
    int map12idx;
    int map13idx;
    map0idx = opDat0Map[n + set_size * 0];
    map10idx = opDat10Map[n + set_size * 0];
    map11idx = opDat10Map[n + set_size * 1];
    map12idx = opDat10Map[n + set_size * 2];
    map13idx = opDat10Map[n + set_size * 3];
    const double* arg10_vec[] = {
       &ind_arg3[map10idx],
       &ind_arg3[map11idx],
       &ind_arg3[map12idx],
       &ind_arg3[map13idx]};
    const double* arg17_vec[] = {
       &ind_arg7[map10idx],
       &ind_arg7[map11idx],
       &ind_arg7[map12idx],
       &ind_arg7[map13idx]};
    const double* arg21_vec[] = {
       &ind_arg8[map10idx],
       &ind_arg8[map11idx],
       &ind_arg8[map12idx],
       &ind_arg8[map13idx]};
    const double* arg25_vec[] = {
       &ind_arg9[map10idx],
       &ind_arg9[map11idx],
       &ind_arg9[map12idx],
       &ind_arg9[map13idx]};

    _pmf_3d_mult_complete_flux_gpu<p, np, npf>(fmask_sh,
                              arg1+n,
                              arg2+n,
                              arg3+n,
                              arg4+n,
                              arg5+n,
                              arg6+n,
                              arg7+n,
                              ind_arg1+map0idx,
                              ind_arg2+map0idx,
                              arg10_vec,
                              ind_arg4+map0idx,
                              ind_arg5+map0idx,
                              ind_arg6+map0idx,
                              arg17_vec,
                              arg21_vec,
                              arg25_vec,
                              ind_arg10+map0idx);
  }
}


//host stub function
void custom_kernel_pmf_3d_mult_complete_flux(const int order, char const *name, op_set set,
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
  op_arg arg14,
  op_arg arg15,
  op_arg arg16,
  op_arg arg17,
  op_arg arg21,
  op_arg arg25,
  op_arg arg29) {

  int nargs = 30;
  op_arg args[30];

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
  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<4; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 20, "double", OP_READ);
  }

  args[14] = arg14;
  args[15] = arg15;
  args[16] = arg16;
  arg17.idx = 0;
  args[17] = arg17;
  for ( int v=1; v<4; v++ ){
    args[17 + v] = op_arg_dat(arg17.dat, v, arg17.map, 20, "double", OP_READ);
  }

  arg21.idx = 0;
  args[21] = arg21;
  for ( int v=1; v<4; v++ ){
    args[21 + v] = op_arg_dat(arg21.dat, v, arg21.map, 20, "double", OP_READ);
  }

  arg25.idx = 0;
  args[25] = arg25;
  for ( int v=1; v<4; v++ ){
    args[25 + v] = op_arg_dat(arg25.dat, v, arg25.map, 20, "double", OP_READ);
  }

  args[29] = arg29;

  int    ninds   = 11;
  int    inds[30] = {0,-1,-1,-1,-1,-1,-1,-1,1,2,3,3,3,3,4,5,6,7,7,7,7,8,8,8,8,9,9,9,9,10};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pmf_3d_mult_complete_flux\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2, 0);
  if (set_size > 0) {

    if (opDat0_custom_pmf_stride_OP2HOST != getSetSizeFromOpArg(&arg0)) {
      opDat0_custom_pmf_stride_OP2HOST = getSetSizeFromOpArg(&arg0);
      cudaMemcpyToSymbol(opDat0_custom_pmf_stride_OP2CONSTANT, &opDat0_custom_pmf_stride_OP2HOST,sizeof(int));
    }
    if (opDat10_custom_pmf_stride_OP2HOST != getSetSizeFromOpArg(&arg10)) {
      opDat10_custom_pmf_stride_OP2HOST = getSetSizeFromOpArg(&arg10);
      cudaMemcpyToSymbol(opDat10_custom_pmf_stride_OP2CONSTANT, &opDat10_custom_pmf_stride_OP2HOST,sizeof(int));
    }
    if (direct_custom_pmf_stride_OP2HOST != getSetSizeFromOpArg(&arg1)) {
      direct_custom_pmf_stride_OP2HOST = getSetSizeFromOpArg(&arg1);
      cudaMemcpyToSymbol(direct_custom_pmf_stride_OP2CONSTANT,&direct_custom_pmf_stride_OP2HOST,sizeof(int));
    }
    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_69
      int nthread = OP_BLOCK_SIZE_69;
    #else
      int nthread = OP_block_size;
    #endif

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_grouped(nargs, args, 2, 0);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        int nblocks = (end-start-1)/nthread+1;
        switch(order) {
          case 1:
            _op_cuda_pmf_3d_mult_complete_flux<1><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg8.data_d,
            (double *)arg9.data_d,
            (double *)arg10.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg17.data_d,
            (double *)arg21.data_d,
            (double *)arg25.data_d,
            (double *)arg29.data_d,
            arg0.map_data_d,
            arg10.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            start,end,set->size+set->exec_size);
          case 2:
            _op_cuda_pmf_3d_mult_complete_flux<2><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg8.data_d,
            (double *)arg9.data_d,
            (double *)arg10.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg17.data_d,
            (double *)arg21.data_d,
            (double *)arg25.data_d,
            (double *)arg29.data_d,
            arg0.map_data_d,
            arg10.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            start,end,set->size+set->exec_size);
          case 3:
            _op_cuda_pmf_3d_mult_complete_flux<3><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg8.data_d,
            (double *)arg9.data_d,
            (double *)arg10.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg17.data_d,
            (double *)arg21.data_d,
            (double *)arg25.data_d,
            (double *)arg29.data_d,
            arg0.map_data_d,
            arg10.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            start,end,set->size+set->exec_size);
        }
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
