template<int p, int dg_npf>
__device__ void _pmf_3d_mult_faces_flux_gpu(const int node, const int *faceNums,
                          const int *fmaskF, const double *nx, const double *ny,
                          const double *nz, const double *fscale, const double *sJ,
                          const double *in, const double **in_p, const double *in_x,
                          const double *in_y, const double *in_z, const double **in_x_p,
                          const double **in_y_p, const double **in_z_p, double *l_x,
                          double *l_y, double *l_z, double *out) {
  if(!(node < DG_NUM_FACES * dg_npf))
    return;
  const int *fmask = &FMASK_cuda[(p - 1) * DG_NUM_FACES * DG_NPF];

  for(int i = 0; i < DG_NUM_FACES; i++) {
    const double gtau = 2.0 * (DG_ORDER + 1) * (DG_ORDER + 2) * fmax(fscale[i * 2], fscale[i * 2 + 1]);
    const double int_fact = 0.5 * sJ[i];
    const double *inR = in_p[i];
    const double *in_xR = in_x_p[i];
    const double *in_yR = in_y_p[i];
    const double *in_zR = in_z_p[i];

    const int findL = faceNums[2 * i] * dg_npf;
    const int *fmaskL = &fmask[faceNums[2 * i] * dg_npf];
    const int *fmaskR = &fmaskF[i * dg_npf];

    const int fmaskIndL = fmaskL[node];
    const int fmaskIndR = fmaskR[node];
    const double diffL_u = in[fmaskIndL] - inR[fmaskIndR];
    const double diffL_u_x = nx[i] * (in_xR[fmaskIndR] + in_x[fmaskIndL]);
    const double diffL_u_y = ny[i] * (in_yR[fmaskIndR] + in_y[fmaskIndL]);
    const double diffL_u_z = nz[i] * (in_zR[fmaskIndR] + in_z[fmaskIndL]);
    const double diffL_u_grad = diffL_u_x + diffL_u_y + diffL_u_z;

    const int indL = findL + node;
    out[indL] = int_fact * (gtau * diffL_u - diffL_u_grad);
    const double l_tmpL = int_fact * -diffL_u;
    l_x[indL] = nx[i] * l_tmpL;
    l_y[indL] = ny[i] * l_tmpL;
    l_z[indL] = nz[i] * l_tmpL;
  }
}

// CUDA kernel function
template<int p>
__global__ void _op_cuda_pmf_3d_mult_faces_flux(
  const int *__restrict ind_arg0,
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  const double *__restrict ind_arg5,
  const double *__restrict ind_arg6,
  const double *__restrict ind_arg7,
  const double *__restrict ind_arg8,
  double *__restrict ind_arg9,
  double *__restrict ind_arg10,
  double *__restrict ind_arg11,
  double *__restrict ind_arg12,
  const int *__restrict opDat0Map,
  const int *__restrict opDat9Map,
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
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int cell = tid / DG_NPF + start;
  int node = tid % DG_NPF;
  if (cell < end) {
    int n = cell;
    //initialise local variables
    int map0idx;
    int map9idx;
    int map10idx;
    int map11idx;
    int map12idx;
    map0idx = opDat0Map[n + set_size * 0];
    map9idx = opDat9Map[n + set_size * 0];
    map10idx = opDat9Map[n + set_size * 1];
    map11idx = opDat9Map[n + set_size * 2];
    map12idx = opDat9Map[n + set_size * 3];
    const double* arg9_vec[] = {
       &ind_arg2[DG_NP * map9idx],
       &ind_arg2[DG_NP * map10idx],
       &ind_arg2[DG_NP * map11idx],
       &ind_arg2[DG_NP * map12idx]};
    const double* arg16_vec[] = {
       &ind_arg6[DG_NP * map9idx],
       &ind_arg6[DG_NP * map10idx],
       &ind_arg6[DG_NP * map11idx],
       &ind_arg6[DG_NP * map12idx]};
    const double* arg20_vec[] = {
       &ind_arg7[DG_NP * map9idx],
       &ind_arg7[DG_NP * map10idx],
       &ind_arg7[DG_NP * map11idx],
       &ind_arg7[DG_NP * map12idx]};
    const double* arg24_vec[] = {
       &ind_arg8[DG_NP * map9idx],
       &ind_arg8[DG_NP * map10idx],
       &ind_arg8[DG_NP * map11idx],
       &ind_arg8[DG_NP * map12idx]};

    // ind_arg0+map0idx*1
    const int npf = (p + 1) * (p + 2) / 2;
    _pmf_3d_mult_faces_flux_gpu<p,npf>(node,
                           arg1 + n * 8,
                           arg2 + n * DG_NUM_FACES * DG_NPF,
                           arg3 + n * 4,
                           arg4 + n * 4,
                           arg5 + n * 4,
                           arg6 + n * 8,
                           arg7 + n * 4,
                           ind_arg1 + map0idx * DG_NP,
                           arg9_vec,
                           ind_arg3 + map0idx * DG_NP,
                           ind_arg4 + map0idx * DG_NP,
                           ind_arg5 + map0idx * DG_NP,
                           arg16_vec,
                           arg20_vec,
                           arg24_vec,
                           ind_arg9 + map0idx * DG_NUM_FACES * DG_NPF,
                           ind_arg10 + map0idx * DG_NUM_FACES * DG_NPF,
                           ind_arg11 + map0idx * DG_NUM_FACES * DG_NPF,
                           ind_arg12 + map0idx * DG_NUM_FACES * DG_NPF);
  }
}


//host stub function
void custom_kernel_pmf_3d_mult_faces_flux(const int order, char const *name, op_set set,
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
  op_arg arg13,
  op_arg arg14,
  op_arg arg15,
  op_arg arg16,
  op_arg arg20,
  op_arg arg24,
  op_arg arg28,
  op_arg arg29,
  op_arg arg30,
  op_arg arg31){

  int nargs = 32;
  op_arg args[32];

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
  for ( int v=1; v<4; v++ ){
    args[9 + v] = op_arg_dat(arg9.dat, v, arg9.map, DG_NP, "double", OP_READ);
  }

  args[13] = arg13;
  args[14] = arg14;
  args[15] = arg15;
  arg16.idx = 0;
  args[16] = arg16;
  for ( int v=1; v<4; v++ ){
    args[16 + v] = op_arg_dat(arg16.dat, v, arg16.map, DG_NP, "double", OP_READ);
  }

  arg20.idx = 0;
  args[20] = arg20;
  for ( int v=1; v<4; v++ ){
    args[20 + v] = op_arg_dat(arg20.dat, v, arg20.map, DG_NP, "double", OP_READ);
  }

  arg24.idx = 0;
  args[24] = arg24;
  for ( int v=1; v<4; v++ ){
    args[24 + v] = op_arg_dat(arg24.dat, v, arg24.map, DG_NP, "double", OP_READ);
  }

  args[28] = arg28;
  args[29] = arg29;
  args[30] = arg30;
  args[31] = arg31;

  int ninds    = 13;
  int inds[32] = {0,-1,-1,-1,-1,-1,-1,-1,1,2,2,2,2,3,4,5,6,6,6,6,7,7,7,7,8,8,8,8,9,10,11,12};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pmf_3d_mult_faces_flux\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {
    //set CUDA execution parameters
    int nthread = OP_block_size;

    for ( int round=0; round<2; round++ ){
      if (round==1) {
        op_mpi_wait_all_grouped(nargs, args, 2);
      }
      int start = round==0 ? 0 : set->core_size;
      int end = round==0 ? set->core_size : set->size + set->exec_size;
      if (end-start>0) {
        int nblocks = ((end-start-1)/nthread+1) * DG_NPF;
        switch(order) {
          case 1:
            _op_cuda_pmf_3d_mult_faces_flux<1><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg8.data_d,
            (double *)arg9.data_d,
            (double *)arg13.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg20.data_d,
            (double *)arg24.data_d,
            (double *)arg28.data_d,
            (double *)arg29.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            arg0.map_data_d,
            arg9.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 2:
            _op_cuda_pmf_3d_mult_faces_flux<2><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg8.data_d,
            (double *)arg9.data_d,
            (double *)arg13.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg20.data_d,
            (double *)arg24.data_d,
            (double *)arg28.data_d,
            (double *)arg29.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            arg0.map_data_d,
            arg9.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 3:
            _op_cuda_pmf_3d_mult_faces_flux<3><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg8.data_d,
            (double *)arg9.data_d,
            (double *)arg13.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg20.data_d,
            (double *)arg24.data_d,
            (double *)arg28.data_d,
            (double *)arg29.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            arg0.map_data_d,
            arg9.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 4:
            _op_cuda_pmf_3d_mult_faces_flux<4><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg8.data_d,
            (double *)arg9.data_d,
            (double *)arg13.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg20.data_d,
            (double *)arg24.data_d,
            (double *)arg28.data_d,
            (double *)arg29.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            arg0.map_data_d,
            arg9.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 5:
            _op_cuda_pmf_3d_mult_faces_flux<5><<<nblocks,nthread>>>(
            (int *)arg0.data_d,
            (double *)arg8.data_d,
            (double *)arg9.data_d,
            (double *)arg13.data_d,
            (double *)arg14.data_d,
            (double *)arg15.data_d,
            (double *)arg16.data_d,
            (double *)arg20.data_d,
            (double *)arg24.data_d,
            (double *)arg28.data_d,
            (double *)arg29.data_d,
            (double *)arg30.data_d,
            (double *)arg31.data_d,
            arg0.map_data_d,
            arg9.map_data_d,
            (int*)arg1.data_d,
            (int*)arg2.data_d,
            (double*)arg3.data_d,
            (double*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            start,end,set->size+set->exec_size);
            break;
        }
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
