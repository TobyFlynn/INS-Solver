#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

#include "dg_compiler_defs.h"

template<int p, int dg_np, int dg_npf>
__device__ void _pmf_3d_mult_cells_part1_gpu(const int ind, const double *mMat, const double *eMat,
                      const double *J, const double *lx, const double *ly, const double *lz,
                      const double *out_tmp, const double *ux, const double *uy, const double *uz,
                      double *outx, double *outy, double *outz, double *out) {
  if(!(ind < dg_np))
    return;
  const double *mmat_mat = &mMat[(p - 1) * 20 * 20];
  const double *emat_mat = &eMat[(p - 1) * 4 * 10 * 20];

  double outx_t = 0.0;
  double outy_t = 0.0;
  double outz_t = 0.0;
  for(int j = 0; j < dg_np; j++) {
    int mat_ind = DG_MAT_IND(ind, j, dg_np, dg_np);
    outx_t += mmat_mat[mat_ind] * ux[j];
    outy_t += mmat_mat[mat_ind] * uy[j];
    outz_t += mmat_mat[mat_ind] * uz[j];
  }
  outx_t *= *J;
  outy_t *= *J;
  outz_t *= *J;
  double out_t = 0.0;
  for(int j = 0; j < 4 * dg_npf; j++) {
    int mat_ind = DG_MAT_IND(ind, j, dg_np, 4 * dg_npf);
    outx_t += emat_mat[mat_ind] * lx[j];
    outy_t += emat_mat[mat_ind] * ly[j];
    outz_t += emat_mat[mat_ind] * lz[j];
    out_t  += emat_mat[mat_ind] * out_tmp[j];
  }
  outx[ind] = outx_t;
  outy[ind] = outy_t;
  outz[ind] = outz_t;
  out[ind] = out_t;
}

template<int p, int dg_np>
__device__ void _pmf_3d_mult_cells_part2_gpu(const int ind, const double *dr,
                            const double *ds, const double *dt,
                            const double *in_r, const double *in_s,
                            const double *in_t, double *out) {
  if(!(ind < dg_np))
    return;
  const double *dr_mat = &dr[(p - 1) * 20 * 20];
  const double *ds_mat = &ds[(p - 1) * 20 * 20];
  const double *dt_mat = &dt[(p - 1) * 20 * 20];

  double tmp = 0.0;
  for(int n = 0; n < dg_np; n++) {
    int mat_ind = DG_MAT_IND(n, ind, dg_np, dg_np);
    tmp += dr_mat[mat_ind] * in_r[n];
    tmp += ds_mat[mat_ind] * in_s[n];
    tmp += dt_mat[mat_ind] * in_t[n];
  }
  out[ind] += tmp;
}

// CUDA kernel function
template<int p, int NUM_CELLS>
__global__ void _op_cuda_pmf_3d_mult_cells_merged(
  const int *__restrict arg0,
  const double *arg1,
  const double *arg2,
  const double *arg3,
  const double *arg4,
  const double *arg5,
  const double *__restrict arg6,
  const double *__restrict arg7,
  const double *__restrict arg8,
  const double *__restrict arg9,
  const double *__restrict arg10,
  const double *__restrict arg11,
  const double *__restrict arg12,
  const double *__restrict arg13,
  const double *__restrict arg14,
  const double *__restrict arg15,
  const double *__restrict arg16,
  const double *__restrict arg17,
  const double *__restrict arg18,
  const double *__restrict arg19,
  const double *__restrict arg20,
  const double *__restrict arg21,
  const double *__restrict arg22,
  double *arg23,
  int   set_size ) {
  __shared__ double ux_shared[NUM_CELLS * 20];
  __shared__ double uy_shared[NUM_CELLS * 20];
  __shared__ double uz_shared[NUM_CELLS * 20];
  __shared__ double lx_shared[NUM_CELLS * 4 * 10];
  __shared__ double ly_shared[NUM_CELLS * 4 * 10];
  __shared__ double lz_shared[NUM_CELLS * 4 * 10];
  __shared__ double out_tmp_shared[NUM_CELLS * 4 * 10];
  __shared__ double tmp_x_shared[NUM_CELLS * 20];
  __shared__ double tmp_y_shared[NUM_CELLS * 20];
  __shared__ double tmp_z_shared[NUM_CELLS * 20];

  //process set elements
  for (int n = threadIdx.x + blockIdx.x * blockDim.x; n - threadIdx.x < set_size * 20; n += blockDim.x * gridDim.x){
    const int node_id = n % 20;
    const int cell_id = n / 20;
    const int local_cell_id = (n / 20) - ((n - threadIdx.x) / 20);
    __syncthreads();
    const int start_ind = ((n - threadIdx.x) / 20) * 20;
    const int num_elem  = min((n - threadIdx.x + blockDim.x) / 20, set_size) - ((n - threadIdx.x) / 20) + 1;
    for(int i = threadIdx.x; i < num_elem * 20; i += blockDim.x) {
      ux_shared[i] = arg20[start_ind + i];
      uy_shared[i] = arg21[start_ind + i];
      uz_shared[i] = arg22[start_ind + i];
    }
    const int start_ind1 = ((n - threadIdx.x) / 20) * 4 * 10;
    for(int i = threadIdx.x; i < num_elem * 4 * 10; i += blockDim.x) {
      lx_shared[i] = arg16[start_ind1 + i];
      ly_shared[i] = arg17[start_ind1 + i];
      lz_shared[i] = arg18[start_ind1 + i];
      out_tmp_shared[i] = arg19[start_ind1 + i];
    }
    __syncthreads();
    //user-supplied kernel call
    const int np  = (p + 1) * (p + 2) * (p + 3) / 6;
    const int npf = (p + 1) * (p + 2) / 2;
    if(n < set_size * 20)
      _pmf_3d_mult_cells_part1_gpu<p,np,npf>(node_id,
                               arg2,
                               arg1,
                               arg15+cell_id*1,
                               lx_shared+local_cell_id*40,
                               ly_shared+local_cell_id*40,
                               lz_shared+local_cell_id*40,
                               out_tmp_shared+local_cell_id*40,
                               ux_shared+local_cell_id*20,
                               uy_shared+local_cell_id*20,
                               uz_shared+local_cell_id*20,
                               tmp_x_shared + local_cell_id * 20,
                               tmp_y_shared + local_cell_id * 20,
                               tmp_z_shared + local_cell_id * 20,
                               arg23+cell_id*20);
    __syncthreads();
    for(int i = threadIdx.x; i < num_elem * 20; i += blockDim.x) {
      int curr_cell = i / 20 + (n - threadIdx.x) / 20;
      ux_shared[i] = *(arg6 + curr_cell) * tmp_x_shared[i] + *(arg9 + curr_cell) * tmp_y_shared[i] + *(arg12 + curr_cell) * tmp_z_shared[i];
      uy_shared[i] = *(arg7 + curr_cell) * tmp_x_shared[i] + *(arg10 + curr_cell) * tmp_y_shared[i] + *(arg13 + curr_cell) * tmp_z_shared[i];
      uz_shared[i] = *(arg8 + curr_cell) * tmp_x_shared[i] + *(arg11 + curr_cell) * tmp_y_shared[i] + *(arg14 + curr_cell) * tmp_z_shared[i];
    }
    __syncthreads();
    if(n < set_size * 20)
      _pmf_3d_mult_cells_part2_gpu<p,np>(node_id, arg3, arg4, arg5,
                                  ux_shared + local_cell_id * 20,
                                  uy_shared + local_cell_id * 20,
                                  uz_shared + local_cell_id * 20,
                                  arg23+cell_id*20);
  }
}


//host stub function
void custom_kernel_pmf_3d_mult_cells_merged(const int order, char const *name, op_set set,
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
  op_arg arg23){

  double*arg1h = (double *)arg1.data;
  double*arg2h = (double *)arg2.data;
  double*arg3h = (double *)arg3.data;
  double*arg4h = (double *)arg4.data;
  double*arg5h = (double *)arg5.data;
  int nargs = 24;
  op_arg args[24];

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

  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pmf_3d_mult_cells_merged");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(2400*sizeof(double));
    consts_bytes += ROUND_UP(1200*sizeof(double));
    consts_bytes += ROUND_UP(1200*sizeof(double));
    consts_bytes += ROUND_UP(1200*sizeof(double));
    consts_bytes += ROUND_UP(1200*sizeof(double));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg1.data   = OP_consts_h + consts_bytes;
    arg1.data_d = OP_consts_d + consts_bytes;
    memcpy(arg1.data, arg1h, 2400*sizeof(double));
    consts_bytes += ROUND_UP(2400*sizeof(double));
    arg2.data   = OP_consts_h + consts_bytes;
    arg2.data_d = OP_consts_d + consts_bytes;
    memcpy(arg2.data, arg2h, 1200*sizeof(double));
    consts_bytes += ROUND_UP(1200*sizeof(double));
    arg3.data   = OP_consts_h + consts_bytes;
    arg3.data_d = OP_consts_d + consts_bytes;
    memcpy(arg3.data, arg3h, 1200*sizeof(double));
    consts_bytes += ROUND_UP(1200*sizeof(double));
    arg4.data   = OP_consts_h + consts_bytes;
    arg4.data_d = OP_consts_d + consts_bytes;
    memcpy(arg4.data, arg4h, 1200*sizeof(double));
    consts_bytes += ROUND_UP(1200*sizeof(double));
    arg5.data   = OP_consts_h + consts_bytes;
    arg5.data_d = OP_consts_d + consts_bytes;
    memcpy(arg5.data, arg5h, 1200*sizeof(double));
    consts_bytes += ROUND_UP(1200*sizeof(double));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    const int nthread = (256 /  20) * 20;
    const int nblocks = 200 < (set->size * 20) / nthread + 1 ? 200 : (set->size * 20) / nthread + 1;
    const int num_cells = (nthread / 20) + 1;

    switch(order) {
      case 1:
        _op_cuda_pmf_3d_mult_cells_merged<1,num_cells><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
          (double *) arg1.data_d,
          (double *) arg2.data_d,
          (double *) arg3.data_d,
          (double *) arg4.data_d,
          (double *) arg5.data_d,
          (double *) arg6.data_d,
          (double *) arg7.data_d,
          (double *) arg8.data_d,
          (double *) arg9.data_d,
          (double *) arg10.data_d,
          (double *) arg11.data_d,
          (double *) arg12.data_d,
          (double *) arg13.data_d,
          (double *) arg14.data_d,
          (double *) arg15.data_d,
          (double *) arg16.data_d,
          (double *) arg17.data_d,
          (double *) arg18.data_d,
          (double *) arg19.data_d,
          (double *) arg20.data_d,
          (double *) arg21.data_d,
          (double *) arg22.data_d,
          (double *) arg23.data_d,
          set->size );
          break;
        case 2:
          _op_cuda_pmf_3d_mult_cells_merged<2,num_cells><<<nblocks,nthread>>>(
            (int *) arg0.data_d,
            (double *) arg1.data_d,
            (double *) arg2.data_d,
            (double *) arg3.data_d,
            (double *) arg4.data_d,
            (double *) arg5.data_d,
            (double *) arg6.data_d,
            (double *) arg7.data_d,
            (double *) arg8.data_d,
            (double *) arg9.data_d,
            (double *) arg10.data_d,
            (double *) arg11.data_d,
            (double *) arg12.data_d,
            (double *) arg13.data_d,
            (double *) arg14.data_d,
            (double *) arg15.data_d,
            (double *) arg16.data_d,
            (double *) arg17.data_d,
            (double *) arg18.data_d,
            (double *) arg19.data_d,
            (double *) arg20.data_d,
            (double *) arg21.data_d,
            (double *) arg22.data_d,
            (double *) arg23.data_d,
            set->size );
            break;
          case 3:
            _op_cuda_pmf_3d_mult_cells_merged<3,num_cells><<<nblocks,nthread>>>(
              (int *) arg0.data_d,
              (double *) arg1.data_d,
              (double *) arg2.data_d,
              (double *) arg3.data_d,
              (double *) arg4.data_d,
              (double *) arg5.data_d,
              (double *) arg6.data_d,
              (double *) arg7.data_d,
              (double *) arg8.data_d,
              (double *) arg9.data_d,
              (double *) arg10.data_d,
              (double *) arg11.data_d,
              (double *) arg12.data_d,
              (double *) arg13.data_d,
              (double *) arg14.data_d,
              (double *) arg15.data_d,
              (double *) arg16.data_d,
              (double *) arg17.data_d,
              (double *) arg18.data_d,
              (double *) arg19.data_d,
              (double *) arg20.data_d,
              (double *) arg21.data_d,
              (double *) arg22.data_d,
              (double *) arg23.data_d,
              set->size );
              break;
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
}
