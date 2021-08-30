//
// auto-generated by op2.py
//

#ifdef _OPENMP
  #include <omp.h>
#endif

// global constants
extern double reynolds;
extern double froude;
extern double weber;
extern double nu0;
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
#include "init_nu_rho_kernel.cpp"
#include "init_cubature_grad_kernel.cpp"
#include "gauss_reverse_kernel.cpp"
#include "init_gauss_grad_kernel.cpp"
#include "init_gauss_grad3_kernel.cpp"
#include "init_gauss_grad4_kernel.cpp"
#include "init_gauss_grad_neighbour_kernel.cpp"
#include "init_gauss_grad5_kernel.cpp"
#include "gauss_gfi_faces2_kernel.cpp"
#include "glb_ind_kernel_kernel.cpp"
#include "glb_ind_kernelBC_kernel.cpp"
#include "poisson_h_kernel.cpp"
#include "poisson_apply_bc_kernel.cpp"
#include "poisson_cells_kernel.cpp"
#include "poisson_edges_kernel.cpp"
#include "poisson_pre_kernel.cpp"
#include "init_gauss_grad3_2_kernel.cpp"
#include "init_gauss_grad4_2_kernel.cpp"
#include "init_gauss_grad5_2_kernel.cpp"
#include "poisson_op1_kernel.cpp"
#include "poisson_op2_kernel.cpp"
#include "poisson_op3_kernel.cpp"
#include "poisson_op4_kernel.cpp"
#include "poisson_op5_kernel.cpp"
#include "pressure_solve_setup_kernel.cpp"
#include "viscosity_solve_setup_kernel.cpp"
#include "set_ic_kernel.cpp"
#include "calc_dt_kernel.cpp"
#include "advection_flux_kernel.cpp"
#include "zero_npf_kernel.cpp"
#include "advection_faces_kernel.cpp"
#include "advection_bc_kernel.cpp"
#include "advection_numerical_flux_kernel.cpp"
#include "advection_surface_tension_kernel.cpp"
#include "advection_intermediate_vel_kernel.cpp"
#include "pressure_mu_kernel.cpp"
#include "zero_g_np1_kernel.cpp"
#include "pressure_bc_kernel.cpp"
#include "pressure_bc2_kernel.cpp"
#include "pressure_rhs_kernel.cpp"
#include "pressure_grad_flux_kernel.cpp"
#include "pressure_update_vel_kernel.cpp"
#include "zero_g_np_kernel.cpp"
#include "viscosity_bc_kernel.cpp"
#include "viscosity_rhs_kernel.cpp"
#include "save_values_kernel.cpp"
#include "calc_h_kernel.cpp"
#include "init_surface_kernel.cpp"
#include "set_rkQ_kernel.cpp"
#include "update_Q_kernel.cpp"
#include "ls_advec_edges_kernel.cpp"
#include "ls_advec_bedges_kernel.cpp"
#include "ls_advec_flux_kernel.cpp"
#include "ls_advec_rhs_kernel.cpp"
#include "ls_sign_kernel.cpp"
#include "ls_flux_kernel.cpp"
#include "ls_bflux_kernel.cpp"
#include "ls_copy_kernel.cpp"
#include "ls_rhs_kernel.cpp"
#include "ls_add_diff_kernel.cpp"
#include "sigma_mult_kernel.cpp"
#include "sigma_flux_kernel.cpp"
#include "sigma_bflux_kernel.cpp"
#include "diff_mult_kernel.cpp"
#include "diff_flux_kernel.cpp"
#include "diff_bflux_kernel.cpp"
#include "ls_reinit_check_kernel.cpp"
#include "ls_step_kernel.cpp"
#include "ls_normalise_kernel.cpp"
#include "ls_local_vis_kernel.cpp"
