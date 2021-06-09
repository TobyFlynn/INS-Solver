//
// auto-generated by op2.py
//

#ifdef _OPENMP
  #include <omp.h>
#endif

// global constants
extern double gam;
extern double mu;
extern double nu;
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
#include "op_lib_cpp.h"

// user kernel files
#include "init_nodes_kernel.cpp"
#include "init_grid_kernel.cpp"
#include "init_edges_kernel.cpp"
#include "init_cubature_grad_kernel.cpp"
#include "init_cubature_kernel.cpp"
#include "init_cubature_OP_kernel.cpp"
#include "gauss_reverse_kernel.cpp"
#include "init_gauss_kernel.cpp"
#include "gauss_tau_kernel.cpp"
#include "gauss_tau_bc_kernel.cpp"
#include "init_gauss_grad_kernel.cpp"
#include "init_gauss_grad2_kernel.cpp"
#include "init_gauss_grad_neighbour_kernel.cpp"
#include "gauss_grad_faces_kernel.cpp"
#include "gauss_op_kernel.cpp"
#include "gauss_gfi_faces_kernel.cpp"
#include "div_kernel.cpp"
#include "curl_kernel.cpp"
#include "grad_kernel.cpp"
#include "cub_grad_kernel.cpp"
#include "cub_grad_weak_kernel.cpp"
#include "cub_div_weak_kernel.cpp"
#include "inv_J_kernel.cpp"
#include "glb_ind_kernel_kernel.cpp"
#include "glb_ind_kernelBC_kernel.cpp"
#include "poisson_mf2_op_kernel.cpp"
#include "poisson_mf2_opf_kernel.cpp"
#include "poisson_mf2_opbf_kernel.cpp"
#include "poisson_mf2_bc_kernel.cpp"
#include "poisson_mf2_apply_bc_kernel.cpp"
#include "poisson_mf2_mass_kernel.cpp"
#include "poisson_mf2_kernel.cpp"
#include "poisson_mf2_faces_kernel.cpp"
#include "poisson_test_init_kernel.cpp"
#include "poisson_test_bc_kernel.cpp"
#include "poisson_test_error_kernel.cpp"
#include "set_ic_kernel.cpp"
#include "calc_dt_kernel.cpp"
#include "advection_flux_kernel.cpp"
#include "advection_faces_kernel.cpp"
#include "advection_bc_kernel.cpp"
#include "advection_numerical_flux_kernel.cpp"
#include "advection_intermediate_vel_kernel.cpp"
#include "pressure_bc_kernel.cpp"
#include "pressure_bc2_kernel.cpp"
#include "pressure_rhs_kernel.cpp"
#include "pressure_update_vel_kernel.cpp"
#include "viscosity_bc_kernel.cpp"
#include "viscosity_rhs_kernel.cpp"
#include "viscosity_reset_bc_kernel.cpp"
#include "lift_drag_kernel.cpp"
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
#include "sigma_flux_kernel.cpp"
#include "sigma_bflux_kernel.cpp"
#include "sigma_mult_kernel.cpp"
#include "diff_flux_kernel.cpp"
#include "diff_bflux_kernel.cpp"
