//
// auto-generated by op2.py
//

//global constants
#ifndef MAX_CONST_SIZE
#define MAX_CONST_SIZE 128
#endif

__constant__ double reynolds_cuda;
__constant__ double froude_cuda;
__constant__ double weber_cuda;
__constant__ double nu0_cuda;
__constant__ double nu1_cuda;
__constant__ double rho0_cuda;
__constant__ double rho1_cuda;
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

//header
#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

void op_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){
  if (!OP_hybrid_gpu) return;
  if (!strcmp(name,"reynolds")) {
    cutilSafeCall(cudaMemcpyToSymbol(reynolds_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"froude")) {
    cutilSafeCall(cudaMemcpyToSymbol(froude_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"weber")) {
    cutilSafeCall(cudaMemcpyToSymbol(weber_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"nu0")) {
    cutilSafeCall(cudaMemcpyToSymbol(nu0_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"nu1")) {
    cutilSafeCall(cudaMemcpyToSymbol(nu1_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"rho0")) {
    cutilSafeCall(cudaMemcpyToSymbol(rho0_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"rho1")) {
    cutilSafeCall(cudaMemcpyToSymbol(rho1_cuda, dat, dim*size));
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
  {
    printf("error: unknown const name\n"); exit(1);
  }
}

//user kernel files
#include "init_nu_rho_kernel.cu"
#include "init_cubature_grad_kernel.cu"
#include "init_cubature_OP_kernel.cu"
#include "gauss_reverse_kernel.cu"
#include "gauss_tau_kernel.cu"
#include "gauss_tau_bc_kernel.cu"
#include "init_gauss_grad_kernel.cu"
#include "init_gauss_grad2_kernel.cu"
#include "init_gauss_grad_neighbour_kernel.cu"
#include "gauss_grad_faces_kernel.cu"
#include "gauss_op_kernel.cu"
#include "gauss_gfi_faces_kernel.cu"
#include "poisson_h_kernel.cu"
#include "poisson_apply_bc_kernel.cu"
#include "poisson_cells_kernel.cu"
#include "poisson_edges_kernel.cu"
#include "poisson_op1_kernel.cu"
#include "poisson_op2_kernel.cu"
#include "poisson_op3_kernel.cu"
#include "poisson_op4_kernel.cu"
#include "poisson_op5_kernel.cu"
#include "pressure_solve_setup_kernel.cu"
#include "viscosity_solve_setup_kernel.cu"
#include "set_ic_kernel.cu"
#include "calc_dt_kernel.cu"
#include "advection_flux_kernel.cu"
#include "advection_faces_kernel.cu"
#include "advection_bc_kernel.cu"
#include "advection_numerical_flux_kernel.cu"
#include "advection_intermediate_vel_kernel.cu"
#include "pressure_mu_kernel.cu"
#include "pressure_bc_kernel.cu"
#include "pressure_bc2_kernel.cu"
#include "pressure_rhs_kernel.cu"
#include "pressure_grad_flux_kernel.cu"
#include "pressure_update_vel_kernel.cu"
#include "viscosity_bc_kernel.cu"
#include "viscosity_rhs_kernel.cu"
#include "viscosity_rhs_rho_kernel.cu"
#include "viscosity_reset_bc_kernel.cu"
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
#include "ls_reinit_check_kernel.cu"
#include "ls_step_kernel.cu"
