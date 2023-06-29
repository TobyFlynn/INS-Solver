template<int p, int dg_npf>
__device__ void _fpmf_3d_mult_faces_flux_gpu(const int node, const int *faceNums,
                          const int *fmaskF, const double *nx, const double *ny,
                          const double *nz, const double *fscale, const double *sJ,
                          const double *tau, const double *in, const double **in_p,
                          const double *factor, const double *in_x,
                          const double *in_y, const double *in_z,
                          const double **in_x_p, const double **in_y_p,
                          const double **in_z_p, double *l_x, double *l_y,
                          double *l_z, double *out) {
  const int *fmask = &FMASK_cuda[(p - 1) * DG_NUM_FACES * DG_NPF];
  if(!(node < DG_NUM_FACES * dg_npf)) return;

  for(int i = 0; i < DG_NUM_FACES; i++) {
    const int findL = faceNums[2 * i] * dg_npf;
    const int *fmaskL = &fmask[faceNums[2 * i] * dg_npf];
    const int *fmaskR = &fmaskF[i * dg_npf];
    const double *inR = in_p[i];
    const double *in_xR = in_x_p[i];
    const double *in_yR = in_y_p[i];
    const double *in_zR = in_z_p[i];
    const double int_fact = 0.5 * sJ[i];

    const int fmaskIndL = fmaskL[node];
    const int fmaskIndR = fmaskR[node];
    const double diffL_u = in[fmaskIndL] - inR[fmaskIndR];
    const double diffL_u_x = nx[i] * (in_xR[fmaskIndR] + in_x[fmaskIndL]);
    const double diffL_u_y = ny[i] * (in_yR[fmaskIndR] + in_y[fmaskIndL]);
    const double diffL_u_z = nz[i] * (in_zR[fmaskIndR] + in_z[fmaskIndL]);
    const double diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

    const int indL = findL + node;
    out[indL] = int_fact * (tau[i] * diffL_u - diffL_u_grad);
    const double l_tmpL = int_fact * factor[fmaskIndL] * -diffL_u;
    l_x[indL] = nx[i] * l_tmpL;
    l_y[indL] = ny[i] * l_tmpL;
    l_z[indL] = nz[i] * l_tmpL;
  }
}

// CUDA kernel function
template<int p>
__global__ void _op_cuda_fpmf_3d_mult_faces_flux(
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
  double *__restrict ind_arg11,
  double *__restrict ind_arg12,
  double *__restrict ind_arg13,
  const int *__restrict opDat0Map,
  const int *__restrict opDat10Map,
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
    map10idx = opDat10Map[n + set_size * 0];
    map11idx = opDat10Map[n + set_size * 1];
    map12idx = opDat10Map[n + set_size * 2];
    map13idx = opDat10Map[n + set_size * 3];
    const double* arg10_vec[] = {
       &ind_arg2[DG_NP * map10idx],
       &ind_arg2[DG_NP * map11idx],
       &ind_arg2[DG_NP * map12idx],
       &ind_arg2[DG_NP * map13idx]};
    const double* arg18_vec[] = {
       &ind_arg7[DG_NP * map10idx],
       &ind_arg7[DG_NP * map11idx],
       &ind_arg7[DG_NP * map12idx],
       &ind_arg7[DG_NP * map13idx]};
    const double* arg22_vec[] = {
       &ind_arg8[DG_NP * map10idx],
       &ind_arg8[DG_NP * map11idx],
       &ind_arg8[DG_NP * map12idx],
       &ind_arg8[DG_NP * map13idx]};
    const double* arg26_vec[] = {
       &ind_arg9[DG_NP * map10idx],
       &ind_arg9[DG_NP * map11idx],
       &ind_arg9[DG_NP * map12idx],
       &ind_arg9[DG_NP * map13idx]};

    // ind_arg0+map0idx*1
    const int npf = (p + 1) * (p + 2) / 2;
    _fpmf_3d_mult_faces_flux_gpu<p,npf>(node,
                            arg1 + n * 8,
                            arg2 + n * DG_NUM_FACES * DG_NPF,
                            arg3 + n * 4,
                            arg4 + n * 4,
                            arg5 + n * 4,
                            arg6 + n * 8,
                            arg7 + n * 4,
                            arg8 + n * 4,
                            ind_arg1 + map0idx * DG_NP,
                            arg10_vec,
                            ind_arg3 + map0idx * DG_NP,
                            ind_arg4 + map0idx * DG_NP,
                            ind_arg5 + map0idx * DG_NP,
                            ind_arg6 + map0idx * DG_NP,
                            arg18_vec,
                            arg22_vec,
                            arg26_vec,
                            ind_arg10 + map0idx * DG_NUM_FACES * DG_NPF,
                            ind_arg11 + map0idx * DG_NUM_FACES * DG_NPF,
                            ind_arg12 + map0idx * DG_NUM_FACES * DG_NPF,
                            ind_arg13 + map0idx * DG_NUM_FACES * DG_NPF);
  }
}


//host stub function
void custom_kernel_fpmf_3d_mult_faces_flux(const int order, char const *name, op_set set,
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
  op_arg arg18,
  op_arg arg22,
  op_arg arg26,
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
  arg10.idx = 0;
  args[10] = arg10;
  for ( int v=1; v<4; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, DG_NP, "double", OP_READ);
  }

  args[14] = arg14;
  args[15] = arg15;
  args[16] = arg16;
  args[17] = arg17;
  arg18.idx = 0;
  args[18] = arg18;
  for ( int v=1; v<4; v++ ){
    args[18 + v] = op_arg_dat(arg18.dat, v, arg18.map, DG_NP, "double", OP_READ);
  }

  arg22.idx = 0;
  args[22] = arg22;
  for ( int v=1; v<4; v++ ){
    args[22 + v] = op_arg_dat(arg22.dat, v, arg22.map, DG_NP, "double", OP_READ);
  }

  arg26.idx = 0;
  args[26] = arg26;
  for ( int v=1; v<4; v++ ){
    args[26 + v] = op_arg_dat(arg26.dat, v, arg26.map, DG_NP, "double", OP_READ);
  }

  args[30] = arg30;
  args[31] = arg31;
  args[32] = arg32;
  args[33] = arg33;

  int ninds    = 14;
  int inds[34] = {0,-1,-1,-1,-1,-1,-1,-1,-1,1,2,2,2,2,3,4,5,6,7,7,7,7,8,8,8,8,9,9,9,9,10,11,12,13};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: fpmf_3d_mult_faces_flux\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2, 0);
  cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_3d_mult_faces_flux<1>, cudaFuncCachePreferL1));
  cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_3d_mult_faces_flux<2>, cudaFuncCachePreferL1));
  cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_3d_mult_faces_flux<3>, cudaFuncCachePreferL1));
  cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_3d_mult_faces_flux<4>, cudaFuncCachePreferL1));
  cutilSafeCall(cudaFuncSetCacheConfig(_op_cuda_fpmf_3d_mult_faces_flux<5>, cudaFuncCachePreferL1));
  if (set_size > 0) {
    //set CUDA execution parameters
    int nthread = OP_block_size;

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_grouped(nargs, args, 2, 0);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        int nblocks = ((end-start-1)/nthread+1) * DG_NPF;
        switch(order) {
          case 1:
            _op_cuda_fpmf_3d_mult_faces_flux<1><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg9.data_d,
            (double *)arg10.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg17.data_d,
            (double *)arg18.data_d,
            (double *)arg22.data_d,
            (double *)arg26.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            (double *)arg32.data_d,
            (double *)arg33.data_d,
            arg0.map_data_d,
            arg10.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            (double*)arg8.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 2:
            _op_cuda_fpmf_3d_mult_faces_flux<2><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg9.data_d,
            (double *)arg10.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg17.data_d,
            (double *)arg18.data_d,
            (double *)arg22.data_d,
            (double *)arg26.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            (double *)arg32.data_d,
            (double *)arg33.data_d,
            arg0.map_data_d,
            arg10.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            (double*)arg8.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 3:
            _op_cuda_fpmf_3d_mult_faces_flux<3><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg9.data_d,
            (double *)arg10.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg17.data_d,
            (double *)arg18.data_d,
            (double *)arg22.data_d,
            (double *)arg26.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            (double *)arg32.data_d,
            (double *)arg33.data_d,
            arg0.map_data_d,
            arg10.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            (double*)arg8.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 4:
            _op_cuda_fpmf_3d_mult_faces_flux<4><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg9.data_d,
            (double *)arg10.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg17.data_d,
            (double *)arg18.data_d,
            (double *)arg22.data_d,
            (double *)arg26.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            (double *)arg32.data_d,
            (double *)arg33.data_d,
            arg0.map_data_d,
            arg10.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            (double*)arg8.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 5:
            _op_cuda_fpmf_3d_mult_faces_flux<5><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg9.data_d,
            (double *)arg10.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg17.data_d,
            (double *)arg18.data_d,
            (double *)arg22.data_d,
            (double *)arg26.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            (double *)arg32.data_d,
            (double *)arg33.data_d,
            arg0.map_data_d,
            arg10.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            (double*)arg8.data_d,
            start,end,set->size+set->exec_size);
            break;
        }
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
