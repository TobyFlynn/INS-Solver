//
// auto-generated by op2.py
//

// global constants
extern double froude;
extern double weber;
extern double nu1;
extern double rho0;
extern double rho1;
extern int FMASK[6];
extern double ic_u;
extern double ic_v;
extern double cubW_g[12];
extern double cubV_g[36];
extern double cubVDr_g[36];
extern double cubVDs_g[36];
extern double gF0Dr_g[9];
extern double gF0Ds_g[9];
extern double gF1Dr_g[9];
extern double gF1Ds_g[9];
extern double gF2Dr_g[9];
extern double gF2Ds_g[9];
extern double gaussW_g[3];
extern double gFInterp0_g[9];
extern double gFInterp1_g[9];
extern double gFInterp2_g[9];
extern double gF0DrR_g[9];
extern double gF0DsR_g[9];
extern double gF1DrR_g[9];
extern double gF1DsR_g[9];
extern double gF2DrR_g[9];
extern double gF2DsR_g[9];
extern double gFInterp0R_g[9];
extern double gFInterp1R_g[9];
extern double gFInterp2R_g[9];

// header
#include "op_lib_c.h"

void op_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){}
// user kernel files
#include "init_nu_rho_acckernel.c"
#include "init_cubature_grad_acckernel.c"
#include "init_cubature_OP_acckernel.c"
#include "gauss_reverse_acckernel.c"
#include "gauss_tau_acckernel.c"
#include "gauss_tau_bc_acckernel.c"
#include "init_gauss_grad_acckernel.c"
#include "init_gauss_grad2_acckernel.c"
#include "init_gauss_grad_neighbour_acckernel.c"
#include "gauss_grad_faces_acckernel.c"
#include "gauss_op_acckernel.c"
#include "gauss_gfi_faces_acckernel.c"
#include "glb_ind_kernel_acckernel.c"
#include "glb_ind_kernelBC_acckernel.c"
#include "poisson_h_acckernel.c"
#include "poisson_apply_bc_acckernel.c"
#include "poisson_cells_acckernel.c"
#include "poisson_edges_acckernel.c"
#include "poisson_pre_acckernel.c"
#include "poisson_op1_acckernel.c"
#include "poisson_op2_acckernel.c"
#include "poisson_op3_acckernel.c"
#include "poisson_op4_acckernel.c"
#include "poisson_op5_acckernel.c"
#include "pressure_solve_setup_acckernel.c"
#include "viscosity_solve_setup_acckernel.c"
#include "set_ic_acckernel.c"
#include "calc_dt_acckernel.c"
#include "advection_flux_acckernel.c"
#include "advection_faces_acckernel.c"
#include "advection_bc_acckernel.c"
#include "advection_numerical_flux_acckernel.c"
#include "advection_intermediate_vel_acckernel.c"
#include "pressure_mu_acckernel.c"
#include "pressure_bc_acckernel.c"
#include "pressure_bc2_acckernel.c"
#include "pressure_rhs_acckernel.c"
#include "pressure_grad_flux_acckernel.c"
#include "pressure_update_vel_acckernel.c"
#include "viscosity_bc_acckernel.c"
#include "viscosity_rhs_acckernel.c"
#include "viscosity_rhs_rho_acckernel.c"
#include "viscosity_reset_bc_acckernel.c"
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
