template<int p, int dg_npf>
__device__ void _fpmf_3d_mult_faces_flux_gpu(const int node, const int *faceNums,
                          const int *fmaskF, const double *nx, const double *ny,
                          const double *nz, const double *fscale, const double *sJ,
                          const double *tau, const double **in, const double **factor,
                          const double **in_x, const double **in_y,
                          const double **in_z, double *l_x, double *l_y,
                          double *l_z, double *out) {
  const int *fmask = &FMASK_cuda[(p - 1) * 4 * 10];
  const double *inL = in[0];
  const double *in_xL = in_x[0];
  const double *in_yL = in_y[0];
  const double *in_zL = in_z[0];
  const double *factorL = factor[0];

  for(int i = 0; i < 4; i++) {
    const int findL = faceNums[2 * i] * dg_npf;
    const int *fmaskL = &fmask[faceNums[2 * i] * dg_npf];
    const int *fmaskR = &fmaskF[i * dg_npf];
    const double *inR = in[i + 1];
    const double *in_xR = in_x[i + 1];
    const double *in_yR = in_y[i + 1];
    const double *in_zR = in_z[i + 1];
    const double *factorR = factor[i + 1];
    const double int_fact = 0.5 * sJ[i];

    // for(int j = 0; j < dg_npf; j++) {
      const int fmaskIndL = fmaskL[node];
      const int fmaskIndR = fmaskR[node];
      const double diffL_u = inL[fmaskIndL] - inR[fmaskIndR];
      const double diffL_u_x = nx[i] * (in_xR[fmaskIndR] + in_xL[fmaskIndL]);
      const double diffL_u_y = ny[i] * (in_yR[fmaskIndR] + in_yL[fmaskIndL]);
      const double diffL_u_z = nz[i] * (in_zR[fmaskIndR] + in_zL[fmaskIndL]);
      const double diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

      const int indL = findL + node;
      out[indL] = int_fact * (tau[i] * diffL_u - diffL_u_grad);
      const double l_tmpL = int_fact * factorL[fmaskIndL] * -diffL_u;
      l_x[indL] = nx[i] * l_tmpL;
      l_y[indL] = ny[i] * l_tmpL;
      l_z[indL] = nz[i] * l_tmpL;
    // }
  }

}

// CUDA kernel function
__global__ void _op_cuda_fpmf_3d_mult_faces_flux(
  const int *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  double *__restrict ind_arg6,
  double *__restrict ind_arg7,
  double *__restrict ind_arg8,
  double *__restrict ind_arg9,
  const int *__restrict opDat0Map,
  const int *__restrict arg1,
  const int *__restrict arg2,
  const double *__restrict arg3,
  const double *__restrict arg4,
  const double *__restrict arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  const double *__restrict arg8,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int cell = tid / DG_NPF + start;
  int node = tid % DG_NPF;
  if (cell < end) {
    int n = cell;
    //initialise local variables
    int map0idx;
    int map10idx;
    int map11idx;
    int map12idx;
    int map13idx;
    map0idx = opDat0Map[n + set_size * 0];
    map10idx = opDat0Map[n + set_size * 1];
    map11idx = opDat0Map[n + set_size * 2];
    map12idx = opDat0Map[n + set_size * 3];
    map13idx = opDat0Map[n + set_size * 4];
    const double* arg9_vec[] = {
       &ind_arg1[20 * map0idx],
       &ind_arg1[20 * map10idx],
       &ind_arg1[20 * map11idx],
       &ind_arg1[20 * map12idx],
       &ind_arg1[20 * map13idx]};
    const double* arg14_vec[] = {
       &ind_arg2[20 * map0idx],
       &ind_arg2[20 * map10idx],
       &ind_arg2[20 * map11idx],
       &ind_arg2[20 * map12idx],
       &ind_arg2[20 * map13idx]};
    const double* arg19_vec[] = {
       &ind_arg3[20 * map0idx],
       &ind_arg3[20 * map10idx],
       &ind_arg3[20 * map11idx],
       &ind_arg3[20 * map12idx],
       &ind_arg3[20 * map13idx]};
    const double* arg24_vec[] = {
       &ind_arg4[20 * map0idx],
       &ind_arg4[20 * map10idx],
       &ind_arg4[20 * map11idx],
       &ind_arg4[20 * map12idx],
       &ind_arg4[20 * map13idx]};
    const double* arg29_vec[] = {
       &ind_arg5[20 * map0idx],
       &ind_arg5[20 * map10idx],
       &ind_arg5[20 * map11idx],
       &ind_arg5[20 * map12idx],
       &ind_arg5[20 * map13idx]};

    //user-supplied kernel call
    // switch(*(ind_arg0+map0idx*1)) {
    //   case 1:
    //     _fpmf_3d_mult_faces_flux_gpu<1,3>(arg1+n*8, arg2+n*40, arg3+n*4, arg4+n*4,
    //               arg5+n*4, arg6+n*8, arg7+n*4, arg8+n*4, arg9_vec, arg14_vec,
    //               arg19_vec, arg24_vec, arg29_vec, ind_arg6+map0idx*40,
    //               ind_arg7+map0idx*40, ind_arg8+map0idx*40, ind_arg9+map0idx*40);
    //     break;
    //   case 2:
    //     _fpmf_3d_mult_faces_flux_gpu<2,6>(arg1+n*8, arg2+n*40, arg3+n*4, arg4+n*4,
    //               arg5+n*4, arg6+n*8, arg7+n*4, arg8+n*4, arg9_vec, arg14_vec,
    //               arg19_vec, arg24_vec, arg29_vec, ind_arg6+map0idx*40,
    //               ind_arg7+map0idx*40, ind_arg8+map0idx*40, ind_arg9+map0idx*40);
    //     break;
    //   case 3:
        _fpmf_3d_mult_faces_flux_gpu<3,10>(node, arg1+n*8, arg2+n*40, arg3+n*4, arg4+n*4,
                  arg5+n*4, arg6+n*8, arg7+n*4, arg8+n*4, arg9_vec, arg14_vec,
                  arg19_vec, arg24_vec, arg29_vec, ind_arg6+map0idx*40,
                  ind_arg7+map0idx*40, ind_arg8+map0idx*40, ind_arg9+map0idx*40);
    //     break;
    // }
  }
}


//host stub function
void custom_kernel_fpmf_3d_mult_faces_flux(char const *name, op_set set,
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
  op_arg arg14,
  op_arg arg19,
  op_arg arg24,
  op_arg arg29,
  op_arg arg34,
  op_arg arg35,
  op_arg arg36,
  op_arg arg37){

  int nargs = 38;
  op_arg args[38];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  arg9.idx = 0;
  args[9] = arg9;
  for ( int v=1; v<5; v++ ){
    args[9 + v] = op_arg_dat(arg9.dat, v, arg9.map, 20, "double", OP_READ);
  }

  arg14.idx = 0;
  args[14] = arg14;
  for ( int v=1; v<5; v++ ){
    args[14 + v] = op_arg_dat(arg14.dat, v, arg14.map, 20, "double", OP_READ);
  }

  arg19.idx = 0;
  args[19] = arg19;
  for ( int v=1; v<5; v++ ){
    args[19 + v] = op_arg_dat(arg19.dat, v, arg19.map, 20, "double", OP_READ);
  }

  arg24.idx = 0;
  args[24] = arg24;
  for ( int v=1; v<5; v++ ){
    args[24 + v] = op_arg_dat(arg24.dat, v, arg24.map, 20, "double", OP_READ);
  }

  arg29.idx = 0;
  args[29] = arg29;
  for ( int v=1; v<5; v++ ){
    args[29 + v] = op_arg_dat(arg29.dat, v, arg29.map, 20, "double", OP_READ);
  }

  args[34] = arg34;
  args[35] = arg35;
  args[36] = arg36;
  args[37] = arg37;

  int    ninds   = 10;
  int    inds[38] = {0,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,7,8,9};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: fpmf_3d_mult_faces_flux\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {
    int nthread = OP_block_size;

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_grouped(nargs, args, 2);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        if(round == 0) {
          int nblocks = ((end-start-1)/nthread+1) * DG_NPF;
          _op_cuda_fpmf_3d_mult_faces_flux<<<nblocks,nthread>>>(
          (int *)arg0.data_d,
          (double *)arg9.data_d,
          (double *)arg14.data_d,
          (double *)arg19.data_d,
          (double *)arg24.data_d,
          (double *)arg29.data_d,
          (double *)arg34.data_d,
          (double *)arg35.data_d,
          (double *)arg36.data_d,
          (double *)arg37.data_d,
          arg0.map_data_d,
          (int*)arg1.data_d,
          (int*)arg2.data_d,
          (double*)arg3.data_d,
          (double*)arg4.data_d,
          (double*)arg5.data_d,
          (double*)arg6.data_d,
          (double*)arg7.data_d,
          (double*)arg8.data_d,
          start,end,set->size+set->exec_size);
        } else {
          int nblocks = (end-start-1)/nthread+1;
          op_cuda_fpmf_3d_mult_faces_flux<<<nblocks,nthread>>>(
          (int *)arg0.data_d,
          (double *)arg9.data_d,
          (double *)arg14.data_d,
          (double *)arg19.data_d,
          (double *)arg24.data_d,
          (double *)arg29.data_d,
          (double *)arg34.data_d,
          (double *)arg35.data_d,
          (double *)arg36.data_d,
          (double *)arg37.data_d,
          arg0.map_data_d,
          (int*)arg1.data_d,
          (int*)arg2.data_d,
          (double*)arg3.data_d,
          (double*)arg4.data_d,
          (double*)arg5.data_d,
          (double*)arg6.data_d,
          (double*)arg7.data_d,
          (double*)arg8.data_d,
          start,end,set->size+set->exec_size);
        }
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
