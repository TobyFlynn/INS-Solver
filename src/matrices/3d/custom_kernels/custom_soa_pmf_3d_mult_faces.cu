__constant__ int opDat0_custom_pmf_stride_OP2CONSTANT;
int opDat0_custom_pmf_stride_OP2HOST=-1;
__constant__ int direct_custom_pmf_stride_OP2CONSTANT;
int direct_custom_pmf_stride_OP2HOST=-1;

template<int p, int npf>
__device__ void _pmf_3d_mult_faces_gpu(const int *faceNum, const int *fmaskL_corrected,
                              const int *fmaskR_corrected, const double *nx,
                              const double *ny, const double *nz, const double *fscale,
                              const double *sJ, const double **in, const double **in_x,
                              const double **in_y, const double **in_z, double **l_x,
                              double **l_y, double **l_z, double **out) {
  const double gtau = 2.0 * (p + 1) * (p + 2) * fmax(fscale[(0)*direct_custom_pmf_stride_OP2CONSTANT], fscale[(1)*direct_custom_pmf_stride_OP2CONSTANT]);

  const double nxL = nx[(0)*direct_custom_pmf_stride_OP2CONSTANT];
  const double nyL = ny[(0)*direct_custom_pmf_stride_OP2CONSTANT];
  const double nzL = nz[(0)*direct_custom_pmf_stride_OP2CONSTANT];
  const double nxR = nx[(1)*direct_custom_pmf_stride_OP2CONSTANT];
  const double nyR = ny[(1)*direct_custom_pmf_stride_OP2CONSTANT];
  const double nzR = nz[(1)*direct_custom_pmf_stride_OP2CONSTANT];
  const double int_factL = 0.5 * sJ[(0)*direct_custom_pmf_stride_OP2CONSTANT];
  const double int_factR = 0.5 * sJ[(0)*direct_custom_pmf_stride_OP2CONSTANT];

  #pragma unroll
  for(int j = 0; j < npf; j++) {
    const int fmaskL_ind = FMASK_cuda[(p - 1) * 4 * 10 + faceNum[(0)*direct_custom_pmf_stride_OP2CONSTANT] * npf + j] * opDat0_custom_pmf_stride_OP2CONSTANT;
    const int fmaskR_ind_corr = fmaskR_corrected[(j)*direct_custom_pmf_stride_OP2CONSTANT] * opDat0_custom_pmf_stride_OP2CONSTANT;
    const double diffL_u = in[0][fmaskL_ind] - in[1][fmaskR_ind_corr];
    const double diffL_u_grad = nxL * (in_x[1][fmaskR_ind_corr] + in_x[0][fmaskL_ind])
                              + nyL * (in_y[1][fmaskR_ind_corr] + in_y[0][fmaskL_ind])
                              + nzL * (in_z[1][fmaskR_ind_corr] + in_z[0][fmaskL_ind]);

    const int indL = (faceNum[(0)*direct_custom_pmf_stride_OP2CONSTANT] * npf + j) * opDat0_custom_pmf_stride_OP2CONSTANT;
    out[0][indL] = int_factL * (gtau * diffL_u - diffL_u_grad);
    const double l_tmpL = int_factL * -diffL_u;
    l_x[0][indL] = nxL * l_tmpL;
    l_y[0][indL] = nyL * l_tmpL;
    l_z[0][indL] = nzL * l_tmpL;

    const int fmaskR_ind = FMASK_cuda[(p - 1) * 4 * 10 + faceNum[(1)*direct_custom_pmf_stride_OP2CONSTANT] * npf + j] * opDat0_custom_pmf_stride_OP2CONSTANT;
    const int fmaskL_ind_corr = fmaskL_corrected[(j)*direct_custom_pmf_stride_OP2CONSTANT] * opDat0_custom_pmf_stride_OP2CONSTANT;
    const double diffR_u = in[1][fmaskR_ind] - in[0][fmaskL_ind_corr];
    const double diffR_u_grad = nxR * (in_x[1][fmaskR_ind] + in_x[0][fmaskL_ind_corr])
                              + nyR * (in_y[1][fmaskR_ind] + in_y[0][fmaskL_ind_corr])
                              + nzR * (in_z[1][fmaskR_ind] + in_z[0][fmaskL_ind_corr]);

    const int indR = (faceNum[(1)*direct_custom_pmf_stride_OP2CONSTANT] * npf + j) * opDat0_custom_pmf_stride_OP2CONSTANT;
    out[1][indR] = int_factR * (gtau * diffR_u - diffR_u_grad);
    const double l_tmpR = int_factR * -diffR_u;
    l_x[1][indR] = nxR * l_tmpR;
    l_y[1][indR] = nyR * l_tmpR;
    l_z[1][indR] = nzR * l_tmpR;
  }

}

template<int p>
__global__ void _op_cuda_pmf_3d_mult_faces(
  const double *__restrict ind_arg1,
  const double *__restrict ind_arg2,
  const double *__restrict ind_arg3,
  const double *__restrict ind_arg4,
  double *__restrict ind_arg5,
  double *__restrict ind_arg6,
  double *__restrict ind_arg7,
  double *__restrict ind_arg8,
  const int *__restrict opDat0Map,
  const int *__restrict arg2,
  const int *__restrict arg3,
  const int *__restrict arg4,
  const double *__restrict arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  const double *__restrict arg8,
  const double *__restrict arg9,
  int start,
  int end,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = tid + start;
    //initialise local variables
    int map0idx;
    int map1idx;
    map0idx = opDat0Map[n + set_size * 0];
    map1idx = opDat0Map[n + set_size * 1];
    const double* arg10_vec[] = {
       &ind_arg1[map0idx],
       &ind_arg1[map1idx]};
    const double* arg12_vec[] = {
       &ind_arg2[map0idx],
       &ind_arg2[map1idx]};
    const double* arg14_vec[] = {
       &ind_arg3[map0idx],
       &ind_arg3[map1idx]};
    const double* arg16_vec[] = {
       &ind_arg4[map0idx],
       &ind_arg4[map1idx]};
    double* arg18_vec[] = {
       &ind_arg5[map0idx],
       &ind_arg5[map1idx]};
    double* arg20_vec[] = {
       &ind_arg6[map0idx],
       &ind_arg6[map1idx]};
    double* arg22_vec[] = {
       &ind_arg7[map0idx],
       &ind_arg7[map1idx]};
    double* arg24_vec[] = {
       &ind_arg8[map0idx],
       &ind_arg8[map1idx]};

    //user-supplied kernel call
    _pmf_3d_mult_faces_gpu<p, (p + 1) * (p + 2) / 2>(
                      arg2+n,
                      arg3+n,
                      arg4+n,
                      arg5+n,
                      arg6+n,
                      arg7+n,
                      arg8+n,
                      arg9+n,
                      arg10_vec,
                      arg12_vec,
                      arg14_vec,
                      arg16_vec,
                      arg18_vec,
                      arg20_vec,
                      arg22_vec,
                      arg24_vec);
  }
}


//host stub function
void custom_kernel_pmf_3d_mult_faces(const int order, char const *name, op_set set,
  op_arg arg0,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9,
  op_arg arg10,
  op_arg arg12,
  op_arg arg14,
  op_arg arg16,
  op_arg arg18,
  op_arg arg20,
  op_arg arg22,
  op_arg arg24) {

  int nargs = 26;
  op_arg args[26];

  arg0.idx = 0;
  args[0] = arg0;
  for ( int v=1; v<2; v++ ){
    args[0 + v] = op_arg_dat(arg0.dat, v, arg0.map, 1, "int", OP_READ);
  }

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
  for ( int v=1; v<2; v++ ){
    args[10 + v] = op_arg_dat(arg10.dat, v, arg10.map, 20, "double", OP_READ);
  }

  arg12.idx = 0;
  args[12] = arg12;
  for ( int v=1; v<2; v++ ){
    args[12 + v] = op_arg_dat(arg12.dat, v, arg12.map, 20, "double", OP_READ);
  }

  arg14.idx = 0;
  args[14] = arg14;
  for ( int v=1; v<2; v++ ){
    args[14 + v] = op_arg_dat(arg14.dat, v, arg14.map, 20, "double", OP_READ);
  }

  arg16.idx = 0;
  args[16] = arg16;
  for ( int v=1; v<2; v++ ){
    args[16 + v] = op_arg_dat(arg16.dat, v, arg16.map, 20, "double", OP_READ);
  }

  arg18.idx = 0;
  args[18] = arg18;
  for ( int v=1; v<2; v++ ){
    args[18 + v] = op_arg_dat(arg18.dat, v, arg18.map, 40, "double", OP_WRITE);
  }

  arg20.idx = 0;
  args[20] = arg20;
  for ( int v=1; v<2; v++ ){
    args[20 + v] = op_arg_dat(arg20.dat, v, arg20.map, 40, "double", OP_WRITE);
  }

  arg22.idx = 0;
  args[22] = arg22;
  for ( int v=1; v<2; v++ ){
    args[22 + v] = op_arg_dat(arg22.dat, v, arg22.map, 40, "double", OP_WRITE);
  }

  arg24.idx = 0;
  args[24] = arg24;
  for ( int v=1; v<2; v++ ){
    args[24 + v] = op_arg_dat(arg24.dat, v, arg24.map, 40, "double", OP_WRITE);
  }

  int    ninds   = 9;
  int    inds[26] = {0,0,-1,-1,-1,-1,-1,-1,-1,-1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: pmf_3d_mult_faces\n");
  }
  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2, 0);
  if (set_size > 0) {

    if ((OP_kernels[69].count==1) || (opDat0_custom_pmf_stride_OP2HOST != getSetSizeFromOpArg(&arg0))) {
      opDat0_custom_pmf_stride_OP2HOST = getSetSizeFromOpArg(&arg0);
      cudaMemcpyToSymbol(opDat0_custom_pmf_stride_OP2CONSTANT, &opDat0_custom_pmf_stride_OP2HOST,sizeof(int));
    }
    if ((OP_kernels[69].count==1) || (direct_custom_pmf_stride_OP2HOST != getSetSizeFromOpArg(&arg2))) {
      direct_custom_pmf_stride_OP2HOST = getSetSizeFromOpArg(&arg2);
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
            _op_cuda_pmf_3d_mult_faces<1><<<nblocks,nthread>>>(
            (double *)arg10.data_d,
            (double *)arg12.data_d,
            (double *)arg14.data_d,
            (double *)arg16.data_d,
            (double *)arg18.data_d,
            (double *)arg20.data_d,
            (double *)arg22.data_d,
            (double *)arg24.data_d,
            arg0.map_data_d,
            (int*)arg2.data_d,
            (int*)arg3.data_d,
            (int*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            (double*)arg8.data_d,
            (double*)arg9.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 2:
            _op_cuda_pmf_3d_mult_faces<2><<<nblocks,nthread>>>(
            (double *)arg10.data_d,
            (double *)arg12.data_d,
            (double *)arg14.data_d,
            (double *)arg16.data_d,
            (double *)arg18.data_d,
            (double *)arg20.data_d,
            (double *)arg22.data_d,
            (double *)arg24.data_d,
            arg0.map_data_d,
            (int*)arg2.data_d,
            (int*)arg3.data_d,
            (int*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            (double*)arg8.data_d,
            (double*)arg9.data_d,
            start,end,set->size+set->exec_size);
            break;
          case 3:
            _op_cuda_pmf_3d_mult_faces<3><<<nblocks,nthread>>>(
            (double *)arg10.data_d,
            (double *)arg12.data_d,
            (double *)arg14.data_d,
            (double *)arg16.data_d,
            (double *)arg18.data_d,
            (double *)arg20.data_d,
            (double *)arg22.data_d,
            (double *)arg24.data_d,
            arg0.map_data_d,
            (int*)arg2.data_d,
            (int*)arg3.data_d,
            (int*)arg4.data_d,
            (double*)arg5.data_d,
            (double*)arg6.data_d,
            (double*)arg7.data_d,
            (double*)arg8.data_d,
            (double*)arg9.data_d,
            start,end,set->size+set->exec_size);
            break;
        }
      }
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
