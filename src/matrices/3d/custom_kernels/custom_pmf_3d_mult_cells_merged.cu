#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

#include "dg_compiler_defs.h"

#include "dg_global_constants/dg_mat_constants_3d.h"

template<int dg_np, int dg_npf>
__device__ void _pmf_3d_mult_cells_part1_gpu(const int ind, const double *mMat, const double *eMat,
                      const double *J, const double *lx, const double *ly, const double *lz,
                      const double *out_tmp, const double *ux, const double *uy, const double *uz,
                      double *outx, double *outy, double *outz, double *out) {
  double outx_t = 0.0;
  double outy_t = 0.0;
  double outz_t = 0.0;
  for(int j = 0; j < dg_np; j++) {
    int mat_ind = DG_MAT_IND(ind, j, dg_np, dg_np);
    outx_t += mMat[mat_ind] * ux[j];
    outy_t += mMat[mat_ind] * uy[j];
    outz_t += mMat[mat_ind] * uz[j];
  }
  outx_t *= *J;
  outy_t *= *J;
  outz_t *= *J;
  double out_t = 0.0;
  for(int j = 0; j < DG_NUM_FACES * dg_npf; j++) {
    int mat_ind = DG_MAT_IND(ind, j, dg_np, DG_NUM_FACES * dg_npf);
    outx_t += eMat[mat_ind] * lx[j];
    outy_t += eMat[mat_ind] * ly[j];
    outz_t += eMat[mat_ind] * lz[j];
    out_t  += eMat[mat_ind] * out_tmp[j];
  }
  outx[ind] = outx_t;
  outy[ind] = outy_t;
  outz[ind] = outz_t;
  out[ind] = out_t;
}

template<int dg_np>
__device__ void _pmf_3d_mult_cells_part2_gpu(const int ind, const double *dr,
                            const double *ds, const double *dt,
                            const double *in_r, const double *in_s,
                            const double *in_t, double *out) {
  double tmp = 0.0;
  for(int n = 0; n < dg_np; n++) {
    int mat_ind = DG_MAT_IND(n, ind, dg_np, dg_np);
    tmp += dr[mat_ind] * in_r[n];
    tmp += ds[mat_ind] * in_s[n];
    tmp += dt[mat_ind] * in_t[n];
  }
  out[ind] += tmp;
}

// CUDA kernel function
template<int p, int NUM_CELLS>
__global__ void _op_cuda_pmf_3d_mult_cells_merged_shared_mat(
  const int *__restrict arg0,
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
  const int np = (p + 1) * (p + 2) * (p + 3) / 6;
  const int npf = (p + 1) * (p + 2) / 2;

  __shared__ double ux_shared[NUM_CELLS * np];
  __shared__ double uy_shared[NUM_CELLS * np];
  __shared__ double uz_shared[NUM_CELLS * np];

  __shared__ double lx_shared[NUM_CELLS * npf * 4];
  __shared__ double ly_shared[NUM_CELLS * npf * 4];
  __shared__ double lz_shared[NUM_CELLS * npf * 4];
  __shared__ double lo_shared[NUM_CELLS * npf * 4];

  __shared__ double tmp_x_shared[NUM_CELLS * np];
  __shared__ double tmp_y_shared[NUM_CELLS * np];
  __shared__ double tmp_z_shared[NUM_CELLS * np];

  __shared__ DG_FP mass_shared[np * np];
  __shared__ DG_FP emat_shared[npf * DG_NUM_FACES * np];
  __shared__ DG_FP dr_shared[np * np];
  __shared__ DG_FP ds_shared[np * np];
  __shared__ DG_FP dt_shared[np * np];

  const int start_ind_mat = (p - 1) * DG_NP * DG_NP;
  for(int i = threadIdx.x; i < np * np; i += blockDim.x) {
    dr_shared[i] = dg_Dr_kernel[start_ind_mat + i];
    ds_shared[i] = dg_Ds_kernel[start_ind_mat + i];
    dt_shared[i] = dg_Dt_kernel[start_ind_mat + i];
    mass_shared[i] = dg_Mass_kernel[start_ind_mat + i];
  }
  const int start_ind_emat = (p - 1) * DG_NP * DG_NPF * DG_NUM_FACES;
  for(int i = threadIdx.x; i < npf * DG_NUM_FACES * np; i += blockDim.x) {
    emat_shared[i] = dg_Emat_kernel[start_ind_emat + i];
  }

  __syncthreads();

  //process set elements
  for (int n = threadIdx.x + blockIdx.x * blockDim.x; n - threadIdx.x < set_size * np; n += blockDim.x * gridDim.x){
    const int node_id = n % np;
    const int cell_id = n / np;
    const int local_cell_id = (n / np) - ((n - threadIdx.x) / np);
    const int start_ind = ((n - threadIdx.x) / np) * DG_NP;
    const int num_elem  = min((n - threadIdx.x + blockDim.x) / np, set_size) - ((n - threadIdx.x) / np) + 1;
    for(int i = threadIdx.x; i < num_elem * np; i += blockDim.x) {
      int curr_cell = i / np + (n - threadIdx.x) / np;
      int curr_node = i % np;
      ux_shared[i] = arg20[curr_cell * DG_NP + curr_node];
      uy_shared[i] = arg21[curr_cell * DG_NP + curr_node];
      uz_shared[i] = arg22[curr_cell * DG_NP + curr_node];
    }
    for(int i = threadIdx.x; i < num_elem * npf * 4; i += blockDim.x) {
      int curr_cell = i / (npf * 4) + (n - threadIdx.x) / np;
      int curr_node = i % (npf * 4);
      lx_shared[i] = arg16[curr_cell * DG_NPF * 4 + curr_node];
      ly_shared[i] = arg17[curr_cell * DG_NPF * 4 + curr_node];
      lz_shared[i] = arg18[curr_cell * DG_NPF * 4 + curr_node];
      lo_shared[i] = arg19[curr_cell * DG_NPF * 4 + curr_node];
    }
    __syncthreads();
    //user-supplied kernel call
    if(n < set_size * np)
      _pmf_3d_mult_cells_part1_gpu<np,npf>(node_id,
                               mass_shared,
                               emat_shared,
                               arg15 + cell_id,
                               lx_shared + local_cell_id * DG_NUM_FACES * npf,
                               ly_shared + local_cell_id * DG_NUM_FACES * npf,
                               lz_shared + local_cell_id * DG_NUM_FACES * npf,
                               lo_shared + local_cell_id * DG_NUM_FACES * npf,
                               ux_shared + local_cell_id * np,
                               uy_shared + local_cell_id * np,
                               uz_shared + local_cell_id * np,
                               tmp_x_shared + local_cell_id * np,
                               tmp_y_shared + local_cell_id * np,
                               tmp_z_shared + local_cell_id * np,
                               arg23 + cell_id * DG_NP);
    __syncthreads();
    for(int i = threadIdx.x; i < num_elem * np; i += blockDim.x) {
      int curr_cell = i / np + (n - threadIdx.x) / np;
      DG_FP tmp_x = tmp_x_shared[i];
      DG_FP tmp_y = tmp_y_shared[i];
      DG_FP tmp_z = tmp_z_shared[i];
      tmp_x_shared[i] = *(arg6 + curr_cell) * tmp_x + *(arg9 + curr_cell) * tmp_y + *(arg12 + curr_cell) * tmp_z;
      tmp_y_shared[i] = *(arg7 + curr_cell) * tmp_x + *(arg10 + curr_cell) * tmp_y + *(arg13 + curr_cell) * tmp_z;
      tmp_z_shared[i] = *(arg8 + curr_cell) * tmp_x + *(arg11 + curr_cell) * tmp_y + *(arg14 + curr_cell) * tmp_z;
    }
    __syncthreads();
    if(n < set_size * np)
      _pmf_3d_mult_cells_part2_gpu<np>(node_id, dr_shared, ds_shared, dt_shared,
                                  tmp_x_shared + local_cell_id * np,
                                  tmp_y_shared + local_cell_id * np,
                                  tmp_z_shared + local_cell_id * np,
                                  arg23 + cell_id * DG_NP);
  }
}

template<int p, int NUM_CELLS>
__global__ void _op_cuda_pmf_3d_mult_cells_merged(
  const int *__restrict arg0,
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
  const int np = (p + 1) * (p + 2) * (p + 3) / 6;
  const int npf = (p + 1) * (p + 2) / 2;

  __shared__ double ux_shared[NUM_CELLS * np];
  __shared__ double uy_shared[NUM_CELLS * np];
  __shared__ double uz_shared[NUM_CELLS * np];

  __shared__ double lx_shared[NUM_CELLS * npf * 4];
  __shared__ double ly_shared[NUM_CELLS * npf * 4];
  __shared__ double lz_shared[NUM_CELLS * npf * 4];
  __shared__ double lo_shared[NUM_CELLS * npf * 4];

  __shared__ double tmp_x_shared[NUM_CELLS * np];
  __shared__ double tmp_y_shared[NUM_CELLS * np];
  __shared__ double tmp_z_shared[NUM_CELLS * np];

  const int start_ind_mat = (p - 1) * DG_NP * DG_NP;
  const int start_ind_emat = (p - 1) * DG_NP * DG_NPF * DG_NUM_FACES;

  __syncthreads();

  //process set elements
  for (int n = threadIdx.x + blockIdx.x * blockDim.x; n - threadIdx.x < set_size * np; n += blockDim.x * gridDim.x){
    const int node_id = n % np;
    const int cell_id = n / np;
    const int local_cell_id = (n / np) - ((n - threadIdx.x) / np);
    const int start_ind = ((n - threadIdx.x) / np) * DG_NP;
    const int num_elem  = min((n - threadIdx.x + blockDim.x) / np, set_size) - ((n - threadIdx.x) / np) + 1;
    for(int i = threadIdx.x; i < num_elem * np; i += blockDim.x) {
      int curr_cell = i / np + (n - threadIdx.x) / np;
      int curr_node = i % np;
      ux_shared[i] = arg20[curr_cell * DG_NP + curr_node];
      uy_shared[i] = arg21[curr_cell * DG_NP + curr_node];
      uz_shared[i] = arg22[curr_cell * DG_NP + curr_node];
    }
    for(int i = threadIdx.x; i < num_elem * npf * 4; i += blockDim.x) {
      int curr_cell = i / (npf * 4) + (n - threadIdx.x) / np;
      int curr_node = i % (npf * 4);
      lx_shared[i] = arg16[curr_cell * DG_NPF * 4 + curr_node];
      ly_shared[i] = arg17[curr_cell * DG_NPF * 4 + curr_node];
      lz_shared[i] = arg18[curr_cell * DG_NPF * 4 + curr_node];
      lo_shared[i] = arg19[curr_cell * DG_NPF * 4 + curr_node];
    }
    __syncthreads();
    //user-supplied kernel call
    if(n < set_size * np)
      _pmf_3d_mult_cells_part1_gpu<np,npf>(node_id,
                               dg_Mass_kernel + start_ind_mat,
                               dg_Emat_kernel + start_ind_emat,
                               arg15 + cell_id,
                               lx_shared + local_cell_id * DG_NUM_FACES * npf,
                               ly_shared + local_cell_id * DG_NUM_FACES * npf,
                               lz_shared + local_cell_id * DG_NUM_FACES * npf,
                               lo_shared + local_cell_id * DG_NUM_FACES * npf,
                               ux_shared + local_cell_id * np,
                               uy_shared + local_cell_id * np,
                               uz_shared + local_cell_id * np,
                               tmp_x_shared + local_cell_id * np,
                               tmp_y_shared + local_cell_id * np,
                               tmp_z_shared + local_cell_id * np,
                               arg23 + cell_id * DG_NP);
    __syncthreads();
    for(int i = threadIdx.x; i < num_elem * np; i += blockDim.x) {
      int curr_cell = i / np + (n - threadIdx.x) / np;
      DG_FP tmp_x = tmp_x_shared[i];
      DG_FP tmp_y = tmp_y_shared[i];
      DG_FP tmp_z = tmp_z_shared[i];
      tmp_x_shared[i] = *(arg6 + curr_cell) * tmp_x + *(arg9 + curr_cell) * tmp_y + *(arg12 + curr_cell) * tmp_z;
      tmp_y_shared[i] = *(arg7 + curr_cell) * tmp_x + *(arg10 + curr_cell) * tmp_y + *(arg13 + curr_cell) * tmp_z;
      tmp_z_shared[i] = *(arg8 + curr_cell) * tmp_x + *(arg11 + curr_cell) * tmp_y + *(arg14 + curr_cell) * tmp_z;
    }
    __syncthreads();
    if(n < set_size * np)
      _pmf_3d_mult_cells_part2_gpu<np>(node_id, dg_Dr_kernel + start_ind_mat,
                  dg_Ds_kernel + start_ind_mat, dg_Dt_kernel + start_ind_mat,
                                  tmp_x_shared + local_cell_id * np,
                                  tmp_y_shared + local_cell_id * np,
                                  tmp_z_shared + local_cell_id * np,
                                  arg23 + cell_id * DG_NP);
  }
}

#include "timing.h"
extern Timing *timer;

//host stub function
void custom_kernel_pmf_3d_mult_cells_merged(const int order, char const *name, op_set set,
  op_arg arg0,
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

  int nargs = 19;
  op_arg args[19];

  args[0] = arg0;
  args[1] = arg6;
  args[2] = arg7;
  args[3] = arg8;
  args[4] = arg9;
  args[5] = arg10;
  args[6] = arg11;
  args[7] = arg12;
  args[8] = arg13;
  args[9] = arg14;
  args[10] = arg15;
  args[11] = arg16;
  args[12] = arg17;
  args[13] = arg18;
  args[14] = arg19;
  args[15] = arg20;
  args[16] = arg21;
  args[17] = arg22;
  args[18] = arg23;

  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  pmf_3d_mult_cells_merged");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {
    //set CUDA execution parameters
    const int np  = (order + 1) * (order + 2) * (order + 3) / 6;
    const int nthread = (256 /  np) * np;
    const int nblocks = 200 < (set->size * np) / nthread + 1 ? 200 : (set->size * np) / nthread + 1;
    // const int num_cells = (nthread / DG_NP) + 1;

    switch(order) {
      case 1:
        _op_cuda_pmf_3d_mult_cells_merged_shared_mat<1,(256 / 4) + 1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
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
        _op_cuda_pmf_3d_mult_cells_merged_shared_mat<2,(256 / 10) + 1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
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
        timer->startTimer("fpmf_cells_merged 3rd order");
        _op_cuda_pmf_3d_mult_cells_merged_shared_mat<3,(256 / 20) + 1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
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
        cutilSafeCall(cudaDeviceSynchronize());
        timer->endTimer("fpmf_cells_merged 3rd order");
        break;
      case 4:
        _op_cuda_pmf_3d_mult_cells_merged<4,(256 / 35) + 1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
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
      case 5:
        _op_cuda_pmf_3d_mult_cells_merged<5,(256 / 56) + 1><<<nblocks,nthread>>>(
          (int *) arg0.data_d,
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
