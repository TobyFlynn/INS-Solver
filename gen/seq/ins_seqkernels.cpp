//
// auto-generated by op2.py
//

// global constants
extern double froude;
extern double weber;
extern double nu1;
extern double rho0;
extern double rho1;
extern int FMASK[12];
extern double ic_u;
extern double ic_v;
extern double cubW_g[36];
extern double cubV_g[360];
extern double cubVDr_g[360];
extern double cubVDs_g[360];
extern double gF0Dr_g[60];
extern double gF0Ds_g[60];
extern double gF1Dr_g[60];
extern double gF1Ds_g[60];
extern double gF2Dr_g[60];
extern double gF2Ds_g[60];
extern double gaussW_g[6];
extern double gFInterp0_g[60];
extern double gFInterp1_g[60];
extern double gFInterp2_g[60];
extern double gF0DrR_g[60];
extern double gF0DsR_g[60];
extern double gF1DrR_g[60];
extern double gF1DsR_g[60];
extern double gF2DrR_g[60];
extern double gF2DsR_g[60];
extern double gFInterp0R_g[60];
extern double gFInterp1R_g[60];
extern double gFInterp2R_g[60];

// header
#include "op_lib_cpp.h"

// user kernel files
#include "init_nu_rho_seqkernel.cpp"
#include "init_cubature_grad_seqkernel.cpp"
#include "gauss_reverse_seqkernel.cpp"
#include "init_gauss_grad_seqkernel.cpp"
#include "init_gauss_grad3_seqkernel.cpp"
#include "init_gauss_grad4_seqkernel.cpp"
#include "init_gauss_grad_neighbour_seqkernel.cpp"
#include "init_gauss_grad5_seqkernel.cpp"
#include "gauss_gfi_faces2_seqkernel.cpp"
#include "glb_ind_kernel_seqkernel.cpp"
#include "glb_ind_kernelBC_seqkernel.cpp"
#include "poisson_h_seqkernel.cpp"
#include "poisson_apply_bc_seqkernel.cpp"
#include "poisson_cells_seqkernel.cpp"
#include "poisson_edges_seqkernel.cpp"
#include "poisson_pre_seqkernel.cpp"
#include "poisson_op1_seqkernel.cpp"
#include "poisson_op2_seqkernel.cpp"
#include "poisson_op3_seqkernel.cpp"
#include "poisson_op4_seqkernel.cpp"
#include "poisson_op5_seqkernel.cpp"
#include "pressure_solve_setup_seqkernel.cpp"
#include "viscosity_solve_setup_seqkernel.cpp"
#include "set_ic_seqkernel.cpp"
#include "calc_dt_seqkernel.cpp"
#include "advection_flux_seqkernel.cpp"
#include "zero_npf_seqkernel.cpp"
#include "advection_faces_seqkernel.cpp"
#include "advection_bc_seqkernel.cpp"
#include "advection_numerical_flux_seqkernel.cpp"
#include "advection_surface_tension_seqkernel.cpp"
#include "advection_intermediate_vel_seqkernel.cpp"
#include "pressure_mu_seqkernel.cpp"
#include "zero_g_np1_seqkernel.cpp"
#include "pressure_bc_seqkernel.cpp"
#include "pressure_bc2_seqkernel.cpp"
#include "pressure_rhs_seqkernel.cpp"
#include "pressure_grad_flux_seqkernel.cpp"
#include "pressure_update_vel_seqkernel.cpp"
#include "zero_g_np_seqkernel.cpp"
#include "viscosity_bc_seqkernel.cpp"
#include "viscosity_rhs_seqkernel.cpp"
#include "save_values_seqkernel.cpp"
#include "calc_h_seqkernel.cpp"
#include "init_surface_seqkernel.cpp"
#include "set_rkQ_seqkernel.cpp"
#include "update_Q_seqkernel.cpp"
#include "ls_advec_edges_seqkernel.cpp"
#include "ls_advec_bedges_seqkernel.cpp"
#include "ls_advec_flux_seqkernel.cpp"
#include "ls_advec_rhs_seqkernel.cpp"
#include "ls_sign_seqkernel.cpp"
#include "ls_flux_seqkernel.cpp"
#include "ls_bflux_seqkernel.cpp"
#include "ls_copy_seqkernel.cpp"
#include "ls_rhs_seqkernel.cpp"
#include "ls_add_diff_seqkernel.cpp"
#include "sigma_mult_seqkernel.cpp"
#include "sigma_flux_seqkernel.cpp"
#include "sigma_bflux_seqkernel.cpp"
#include "diff_mult_seqkernel.cpp"
#include "diff_flux_seqkernel.cpp"
#include "diff_bflux_seqkernel.cpp"
#include "ls_reinit_check_seqkernel.cpp"
#include "ls_step_seqkernel.cpp"
#include "ls_normalise_seqkernel.cpp"
#include "ls_local_vis_seqkernel.cpp"
