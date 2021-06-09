//
// auto-generated by op2.py
//

// global constants
double gam_ompkernel;
double mu_ompkernel;
double nu0_ompkernel;
double nu1_ompkernel;
double bc_mach_ompkernel;
double bc_alpha_ompkernel;
double bc_p_ompkernel;
double bc_u_ompkernel;
double bc_v_ompkernel;
int FMASK_ompkernel[15];
double ic_u_ompkernel;
double ic_v_ompkernel;
double cubW_g_ompkernel[46];
double cubV_g_ompkernel[690];
double cubVDr_g_ompkernel[690];
double cubVDs_g_ompkernel[690];
double gF0Dr_g_ompkernel[105];
double gF0Ds_g_ompkernel[105];
double gF1Dr_g_ompkernel[105];
double gF1Ds_g_ompkernel[105];
double gF2Dr_g_ompkernel[105];
double gF2Ds_g_ompkernel[105];
double gaussW_g_ompkernel[7];
double gFInterp0_g_ompkernel[105];
double gFInterp1_g_ompkernel[105];
double gFInterp2_g_ompkernel[105];
double gF0DrR_g_ompkernel[105];
double gF0DsR_g_ompkernel[105];
double gF1DrR_g_ompkernel[105];
double gF1DsR_g_ompkernel[105];
double gF2DrR_g_ompkernel[105];
double gF2DsR_g_ompkernel[105];
double gFInterp0R_g_ompkernel[105];
double gFInterp1R_g_ompkernel[105];
double gFInterp2R_g_ompkernel[105];
double lift_drag_vec_ompkernel[5];

// header
#include "op_lib_cpp.h"

void op_decl_const_char(int dim, char const *type,
  int size, char *dat, char const *name){
  if(!strcmp(name, "gam")) {
    memcpy(&gam_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gam_ompkernel)
  } else if(!strcmp(name, "mu")) {
    memcpy(&mu_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:mu_ompkernel)
  } else if(!strcmp(name, "nu0")) {
    memcpy(&nu0_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:nu0_ompkernel)
  } else if(!strcmp(name, "nu1")) {
    memcpy(&nu1_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:nu1_ompkernel)
  } else if(!strcmp(name, "bc_mach")) {
    memcpy(&bc_mach_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_mach_ompkernel)
  } else if(!strcmp(name, "bc_alpha")) {
    memcpy(&bc_alpha_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_alpha_ompkernel)
  } else if(!strcmp(name, "bc_p")) {
    memcpy(&bc_p_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_p_ompkernel)
  } else if(!strcmp(name, "bc_u")) {
    memcpy(&bc_u_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_u_ompkernel)
  } else if(!strcmp(name, "bc_v")) {
    memcpy(&bc_v_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_v_ompkernel)
  } else if(!strcmp(name, "FMASK")) {
    memcpy(FMASK_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:FMASK_ompkernel[:15])
  } else if(!strcmp(name, "ic_u")) {
    memcpy(&ic_u_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:ic_u_ompkernel)
  } else if(!strcmp(name, "ic_v")) {
    memcpy(&ic_v_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:ic_v_ompkernel)
  } else if(!strcmp(name, "cubW_g")) {
    memcpy(cubW_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:cubW_g_ompkernel[:46])
  } else if(!strcmp(name, "cubV_g")) {
    memcpy(cubV_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:cubV_g_ompkernel[:690])
  } else if(!strcmp(name, "cubVDr_g")) {
    memcpy(cubVDr_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:cubVDr_g_ompkernel[:690])
  } else if(!strcmp(name, "cubVDs_g")) {
    memcpy(cubVDs_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:cubVDs_g_ompkernel[:690])
  } else if(!strcmp(name, "gF0Dr_g")) {
    memcpy(gF0Dr_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF0Dr_g_ompkernel[:105])
  } else if(!strcmp(name, "gF0Ds_g")) {
    memcpy(gF0Ds_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF0Ds_g_ompkernel[:105])
  } else if(!strcmp(name, "gF1Dr_g")) {
    memcpy(gF1Dr_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF1Dr_g_ompkernel[:105])
  } else if(!strcmp(name, "gF1Ds_g")) {
    memcpy(gF1Ds_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF1Ds_g_ompkernel[:105])
  } else if(!strcmp(name, "gF2Dr_g")) {
    memcpy(gF2Dr_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF2Dr_g_ompkernel[:105])
  } else if(!strcmp(name, "gF2Ds_g")) {
    memcpy(gF2Ds_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF2Ds_g_ompkernel[:105])
  } else if(!strcmp(name, "gaussW_g")) {
    memcpy(gaussW_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gaussW_g_ompkernel[:7])
  } else if(!strcmp(name, "gFInterp0_g")) {
    memcpy(gFInterp0_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gFInterp0_g_ompkernel[:105])
  } else if(!strcmp(name, "gFInterp1_g")) {
    memcpy(gFInterp1_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gFInterp1_g_ompkernel[:105])
  } else if(!strcmp(name, "gFInterp2_g")) {
    memcpy(gFInterp2_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gFInterp2_g_ompkernel[:105])
  } else if(!strcmp(name, "gF0DrR_g")) {
    memcpy(gF0DrR_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF0DrR_g_ompkernel[:105])
  } else if(!strcmp(name, "gF0DsR_g")) {
    memcpy(gF0DsR_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF0DsR_g_ompkernel[:105])
  } else if(!strcmp(name, "gF1DrR_g")) {
    memcpy(gF1DrR_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF1DrR_g_ompkernel[:105])
  } else if(!strcmp(name, "gF1DsR_g")) {
    memcpy(gF1DsR_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF1DsR_g_ompkernel[:105])
  } else if(!strcmp(name, "gF2DrR_g")) {
    memcpy(gF2DrR_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF2DrR_g_ompkernel[:105])
  } else if(!strcmp(name, "gF2DsR_g")) {
    memcpy(gF2DsR_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gF2DsR_g_ompkernel[:105])
  } else if(!strcmp(name, "gFInterp0R_g")) {
    memcpy(gFInterp0R_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gFInterp0R_g_ompkernel[:105])
  } else if(!strcmp(name, "gFInterp1R_g")) {
    memcpy(gFInterp1R_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gFInterp1R_g_ompkernel[:105])
  } else if(!strcmp(name, "gFInterp2R_g")) {
    memcpy(gFInterp2R_g_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gFInterp2R_g_ompkernel[:105])
  } else if(!strcmp(name, "lift_drag_vec")) {
    memcpy(lift_drag_vec_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:lift_drag_vec_ompkernel[:5])
  }
}
// user kernel files
#include "init_nodes_omp4kernel_func.cpp"
#include "init_grid_omp4kernel_func.cpp"
#include "init_edges_omp4kernel_func.cpp"
#include "init_cubature_grad_omp4kernel_func.cpp"
#include "init_cubature_omp4kernel_func.cpp"
#include "init_cubature_OP_omp4kernel_func.cpp"
#include "gauss_reverse_omp4kernel_func.cpp"
#include "init_gauss_omp4kernel_func.cpp"
#include "gauss_tau_omp4kernel_func.cpp"
#include "gauss_tau_bc_omp4kernel_func.cpp"
#include "init_gauss_grad_omp4kernel_func.cpp"
#include "init_gauss_grad2_omp4kernel_func.cpp"
#include "init_gauss_grad_neighbour_omp4kernel_func.cpp"
#include "gauss_grad_faces_omp4kernel_func.cpp"
#include "gauss_op_omp4kernel_func.cpp"
#include "gauss_gfi_faces_omp4kernel_func.cpp"
#include "div_omp4kernel_func.cpp"
#include "curl_omp4kernel_func.cpp"
#include "grad_omp4kernel_func.cpp"
#include "cub_grad_omp4kernel_func.cpp"
#include "cub_grad_weak_omp4kernel_func.cpp"
#include "cub_div_weak_omp4kernel_func.cpp"
#include "inv_J_omp4kernel_func.cpp"
#include "glb_ind_kernel_omp4kernel_func.cpp"
#include "glb_ind_kernelBC_omp4kernel_func.cpp"
#include "poisson_mf2_op_omp4kernel_func.cpp"
#include "poisson_mf2_opf_omp4kernel_func.cpp"
#include "poisson_mf2_opbf_omp4kernel_func.cpp"
#include "poisson_mf2_bc_omp4kernel_func.cpp"
#include "poisson_mf2_apply_bc_omp4kernel_func.cpp"
#include "poisson_mf2_mass_omp4kernel_func.cpp"
#include "poisson_mf2_omp4kernel_func.cpp"
#include "poisson_mf2_faces_omp4kernel_func.cpp"
#include "poisson_test_init_omp4kernel_func.cpp"
#include "poisson_test_bc_omp4kernel_func.cpp"
#include "poisson_test_error_omp4kernel_func.cpp"
#include "set_ic_omp4kernel_func.cpp"
#include "calc_dt_omp4kernel_func.cpp"
#include "advection_flux_omp4kernel_func.cpp"
#include "advection_faces_omp4kernel_func.cpp"
#include "advection_bc_omp4kernel_func.cpp"
#include "advection_numerical_flux_omp4kernel_func.cpp"
#include "advection_intermediate_vel_omp4kernel_func.cpp"
#include "pressure_bc_omp4kernel_func.cpp"
#include "pressure_bc2_omp4kernel_func.cpp"
#include "pressure_rhs_omp4kernel_func.cpp"
#include "pressure_update_vel_omp4kernel_func.cpp"
#include "viscosity_bc_omp4kernel_func.cpp"
#include "viscosity_rhs_omp4kernel_func.cpp"
#include "viscosity_reset_bc_omp4kernel_func.cpp"
#include "lift_drag_omp4kernel_func.cpp"
#include "save_values_omp4kernel_func.cpp"
#include "calc_h_omp4kernel_func.cpp"
#include "init_surface_omp4kernel_func.cpp"
#include "set_rkQ_omp4kernel_func.cpp"
#include "update_Q_omp4kernel_func.cpp"
#include "ls_advec_edges_omp4kernel_func.cpp"
#include "ls_advec_bedges_omp4kernel_func.cpp"
#include "ls_advec_flux_omp4kernel_func.cpp"
#include "ls_advec_rhs_omp4kernel_func.cpp"
#include "ls_sign_omp4kernel_func.cpp"
#include "ls_flux_omp4kernel_func.cpp"
#include "ls_bflux_omp4kernel_func.cpp"
#include "ls_copy_omp4kernel_func.cpp"
#include "ls_rhs_omp4kernel_func.cpp"
#include "ls_add_diff_omp4kernel_func.cpp"
#include "sigma_flux_omp4kernel_func.cpp"
#include "sigma_bflux_omp4kernel_func.cpp"
#include "sigma_mult_omp4kernel_func.cpp"
#include "diff_flux_omp4kernel_func.cpp"
#include "diff_bflux_omp4kernel_func.cpp"
#include "ls_reinit_check_omp4kernel_func.cpp"
#include "ls_step_omp4kernel_func.cpp"
