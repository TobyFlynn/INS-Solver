//
// auto-generated by op2.py
//

//global constants
#ifndef MAX_CONST_SIZE
#define MAX_CONST_SIZE 128
#endif

__constant__ double gam_cuda;
__constant__ double mu_cuda;
__constant__ double nu_cuda;
__constant__ double bc_mach_cuda;
__constant__ double bc_alpha_cuda;
__constant__ double bc_p_cuda;
__constant__ double bc_u_cuda;
__constant__ double bc_v_cuda;
__constant__ int FMASK_cuda[15];
__constant__ double ic_u_cuda;
__constant__ double ic_v_cuda;
__constant__ double cubW_g_cuda[46];
__constant__ double cubV_g_cuda[690];
__constant__ double cubVDr_g_cuda[690];
__constant__ double cubVDs_g_cuda[690];
__constant__ double gF0Dr_g_cuda[105];
__constant__ double gF0Ds_g_cuda[105];
__constant__ double gF1Dr_g_cuda[105];
__constant__ double gF1Ds_g_cuda[105];
__constant__ double gF2Dr_g_cuda[105];
__constant__ double gF2Ds_g_cuda[105];
__constant__ double gaussW_g_cuda[7];
__constant__ double gFInterp0_g_cuda[105];
__constant__ double gFInterp1_g_cuda[105];
__constant__ double gFInterp2_g_cuda[105];
__constant__ double gF0DrR_g_cuda[105];
__constant__ double gF0DsR_g_cuda[105];
__constant__ double gF1DrR_g_cuda[105];
__constant__ double gF1DsR_g_cuda[105];
__constant__ double gF2DrR_g_cuda[105];
__constant__ double gF2DsR_g_cuda[105];
__constant__ double gFInterp0R_g_cuda[105];
__constant__ double gFInterp1R_g_cuda[105];
__constant__ double gFInterp2R_g_cuda[105];
__constant__ double lift_drag_vec_cuda[5];

//header
#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

void op_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){
  if (!OP_hybrid_gpu) return;
  if (!strcmp(name,"gam")) {
    cutilSafeCall(cudaMemcpyToSymbol(gam_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"mu")) {
    cutilSafeCall(cudaMemcpyToSymbol(mu_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"nu")) {
    cutilSafeCall(cudaMemcpyToSymbol(nu_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_mach")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_mach_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_alpha")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_alpha_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_p")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_p_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_u")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_u_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_v")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_v_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"FMASK")) {
    cutilSafeCall(cudaMemcpyToSymbol(FMASK_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"ic_u")) {
    cutilSafeCall(cudaMemcpyToSymbol(ic_u_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"ic_v")) {
    cutilSafeCall(cudaMemcpyToSymbol(ic_v_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"cubW_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(cubW_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"cubV_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(cubV_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"cubVDr_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(cubVDr_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"cubVDs_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(cubVDs_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF0Dr_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF0Dr_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF0Ds_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF0Ds_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF1Dr_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF1Dr_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF1Ds_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF1Ds_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF2Dr_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF2Dr_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF2Ds_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF2Ds_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gaussW_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gaussW_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gFInterp0_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gFInterp0_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gFInterp1_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gFInterp1_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gFInterp2_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gFInterp2_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF0DrR_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF0DrR_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF0DsR_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF0DsR_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF1DrR_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF1DrR_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF1DsR_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF1DsR_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF2DrR_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF2DrR_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gF2DsR_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gF2DsR_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gFInterp0R_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gFInterp0R_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gFInterp1R_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gFInterp1R_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"gFInterp2R_g")) {
    cutilSafeCall(cudaMemcpyToSymbol(gFInterp2R_g_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"lift_drag_vec")) {
    cutilSafeCall(cudaMemcpyToSymbol(lift_drag_vec_cuda, dat, dim*size));
  }
  else
  {
    printf("error: unknown const name\n"); exit(1);
  }
}

//user kernel files
#include "init_nodes_kernel.cu"
#include "init_grid_kernel.cu"
#include "init_cubature_grad_kernel.cu"
#include "init_cubature_kernel.cu"
#include "init_cubature_OP_kernel.cu"
#include "gauss_reverse_kernel.cu"
#include "init_gauss_kernel.cu"
#include "gauss_tau_kernel.cu"
#include "gauss_tau_bc_kernel.cu"
#include "init_gauss_grad_kernel.cu"
#include "init_gauss_grad2_kernel.cu"
#include "init_gauss_grad_neighbour_kernel.cu"
#include "gauss_grad_faces_kernel.cu"
#include "gauss_op_kernel.cu"
#include "gauss_gfi_faces_kernel.cu"
#include "div_kernel.cu"
#include "curl_kernel.cu"
#include "grad_kernel.cu"
#include "cub_grad_kernel.cu"
#include "cub_grad_weak_kernel.cu"
#include "cub_div_weak_kernel.cu"
#include "inv_J_kernel.cu"
#include "glb_ind_kernel_kernel.cu"
#include "glb_ind_kernelBC_kernel.cu"
#include "poisson_mf2_op_kernel.cu"
#include "poisson_mf2_opf_kernel.cu"
#include "poisson_mf2_opbf_kernel.cu"
#include "poisson_mf2_bc_kernel.cu"
#include "poisson_mf2_apply_bc_kernel.cu"
#include "poisson_mf2_mass_kernel.cu"
#include "poisson_mf2_kernel.cu"
#include "poisson_mf2_faces_kernel.cu"
#include "poisson_test_init_kernel.cu"
#include "poisson_test_bc_kernel.cu"
#include "poisson_test_error_kernel.cu"
#include "set_ic_kernel.cu"
#include "calc_dt_kernel.cu"
#include "advection_flux_kernel.cu"
#include "advection_faces_kernel.cu"
#include "advection_bc_kernel.cu"
#include "advection_numerical_flux_kernel.cu"
#include "advection_intermediate_vel_kernel.cu"
#include "pressure_bc_kernel.cu"
#include "pressure_bc2_kernel.cu"
#include "pressure_rhs_kernel.cu"
#include "pressure_update_vel_kernel.cu"
#include "viscosity_bc_kernel.cu"
#include "viscosity_rhs_kernel.cu"
#include "viscosity_reset_bc_kernel.cu"
#include "lift_drag_kernel.cu"
#include "save_values_kernel.cu"
#include "calc_h_kernel.cu"
#include "init_surface_kernel.cu"
#include "set_rkQ_kernel.cu"
#include "update_Q_kernel.cu"
#include "ls_advec_edges_kernel.cu"
#include "ls_advec_bedges_kernel.cu"
#include "ls_advec_flux_kernel.cu"
#include "ls_advec_rhs_kernel.cu"
#include "ls_sign_kernel.cu"
#include "ls_flux_kernel.cu"
#include "ls_bflux_kernel.cu"
#include "ls_copy_kernel.cu"
#include "ls_rhs_kernel.cu"
#include "ls_add_diff_kernel.cu"
#include "sigma_flux_kernel.cu"
#include "sigma_bflux_kernel.cu"
#include "sigma_mult_kernel.cu"
#include "diff_flux_kernel.cu"
#include "diff_bflux_kernel.cu"
