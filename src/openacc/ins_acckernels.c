//
// auto-generated by op2.py
//

// global constants
extern double gam;
extern double mu;
extern double nu0;
extern double nu1;
extern double bc_mach;
extern double bc_alpha;
extern double bc_p;
extern double bc_u;
extern double bc_v;
extern int FMASK[15];
extern double ic_u;
extern double ic_v;
extern double cubW_g[46];
extern double cubV_g[690];
extern double cubVDr_g[690];
extern double cubVDs_g[690];
extern double gF0Dr_g[105];
extern double gF0Ds_g[105];
extern double gF1Dr_g[105];
extern double gF1Ds_g[105];
extern double gF2Dr_g[105];
extern double gF2Ds_g[105];
extern double gaussW_g[7];
extern double gFInterp0_g[105];
extern double gFInterp1_g[105];
extern double gFInterp2_g[105];
extern double gF0DrR_g[105];
extern double gF0DsR_g[105];
extern double gF1DrR_g[105];
extern double gF1DsR_g[105];
extern double gF2DrR_g[105];
extern double gF2DsR_g[105];
extern double gFInterp0R_g[105];
extern double gFInterp1R_g[105];
extern double gFInterp2R_g[105];
extern double lift_drag_vec[5];

// header
#include "op_lib_c.h"

void op_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){}
// user kernel files
#include "init_nodes_acckernel.c"
#include "init_grid_acckernel.c"
#include "init_edges_acckernel.c"
#include "init_cubature_grad_acckernel.c"
#include "init_cubature_acckernel.c"
#include "init_cubature_OP_acckernel.c"
#include "gauss_reverse_acckernel.c"
#include "init_gauss_acckernel.c"
#include "gauss_tau_acckernel.c"
#include "gauss_tau_bc_acckernel.c"
#include "init_gauss_grad_acckernel.c"
#include "init_gauss_grad2_acckernel.c"
#include "init_gauss_grad_neighbour_acckernel.c"
#include "gauss_grad_faces_acckernel.c"
#include "gauss_op_acckernel.c"
#include "gauss_gfi_faces_acckernel.c"
#include "div_acckernel.c"
#include "curl_acckernel.c"
#include "grad_acckernel.c"
#include "cub_grad_acckernel.c"
#include "cub_grad_weak_acckernel.c"
#include "cub_div_weak_acckernel.c"
#include "inv_J_acckernel.c"
#include "glb_ind_kernel_acckernel.c"
#include "glb_ind_kernelBC_acckernel.c"
#include "poisson_mf2_op_acckernel.c"
#include "poisson_mf2_opf_acckernel.c"
#include "poisson_mf2_opbf_acckernel.c"
#include "poisson_mf2_bc_acckernel.c"
#include "poisson_mf2_apply_bc_vis_acckernel.c"
#include "poisson_mf2_apply_bc_acckernel.c"
#include "poisson_mf2_mass_acckernel.c"
#include "poisson_mf2_vis_acckernel.c"
#include "poisson_mf2_faces_vis_acckernel.c"
#include "poisson_mf2_acckernel.c"
#include "poisson_mf2_faces_acckernel.c"
#include "poisson_test_init_acckernel.c"
#include "poisson_test_bc_acckernel.c"
#include "poisson_test_error_acckernel.c"
#include "set_ic_acckernel.c"
#include "calc_dt_acckernel.c"
#include "advection_flux_acckernel.c"
#include "advection_faces_acckernel.c"
#include "advection_bc_acckernel.c"
#include "advection_numerical_flux_acckernel.c"
#include "advection_intermediate_vel_acckernel.c"
#include "pressure_bc_acckernel.c"
#include "pressure_bc2_acckernel.c"
#include "pressure_rhs_acckernel.c"
#include "pressure_update_vel_acckernel.c"
#include "viscosity_bc_acckernel.c"
#include "viscosity_rhs_acckernel.c"
#include "viscosity_reset_bc_acckernel.c"
#include "lift_drag_acckernel.c"
#include "save_values_acckernel.c"
#include "calc_h_acckernel.c"
#include "init_surface_acckernel.c"
#include "set_rkQ_acckernel.c"
#include "update_Q_acckernel.c"
#include "ls_advec_edges_acckernel.c"
#include "ls_advec_bedges_acckernel.c"
#include "ls_advec_flux_acckernel.c"
#include "ls_advec_rhs_acckernel.c"
#include "ls_sign_acckernel.c"
#include "ls_flux_acckernel.c"
#include "ls_bflux_acckernel.c"
#include "ls_copy_acckernel.c"
#include "ls_rhs_acckernel.c"
#include "ls_add_diff_acckernel.c"
#include "sigma_flux_acckernel.c"
#include "sigma_bflux_acckernel.c"
#include "sigma_mult_acckernel.c"
#include "diff_flux_acckernel.c"
#include "diff_bflux_acckernel.c"
#include "ls_reinit_check_acckernel.c"
#include "ls_step_acckernel.c"
