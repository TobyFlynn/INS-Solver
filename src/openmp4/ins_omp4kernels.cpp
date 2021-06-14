//
// auto-generated by op2.py
//

// header
#include "op_lib_cpp.h"

// user kernel files
#include "init_nodes_omp4kernel.cpp"
#include "init_grid_omp4kernel.cpp"
#include "init_edges_omp4kernel.cpp"
#include "init_cubature_grad_omp4kernel.cpp"
#include "init_cubature_omp4kernel.cpp"
#include "init_cubature_OP_omp4kernel.cpp"
#include "gauss_reverse_omp4kernel.cpp"
#include "init_gauss_omp4kernel.cpp"
#include "gauss_tau_omp4kernel.cpp"
#include "gauss_tau_bc_omp4kernel.cpp"
#include "init_gauss_grad_omp4kernel.cpp"
#include "init_gauss_grad2_omp4kernel.cpp"
#include "init_gauss_grad_neighbour_omp4kernel.cpp"
#include "gauss_grad_faces_omp4kernel.cpp"
#include "gauss_op_omp4kernel.cpp"
#include "gauss_gfi_faces_omp4kernel.cpp"
#include "div_omp4kernel.cpp"
#include "curl_omp4kernel.cpp"
#include "grad_omp4kernel.cpp"
#include "cub_grad_omp4kernel.cpp"
#include "cub_div_omp4kernel.cpp"
#include "cub_grad_weak_omp4kernel.cpp"
#include "cub_div_weak_omp4kernel.cpp"
#include "inv_J_omp4kernel.cpp"
#include "glb_ind_kernel_omp4kernel.cpp"
#include "glb_ind_kernelBC_omp4kernel.cpp"
#include "poisson_mf2_op_omp4kernel.cpp"
#include "poisson_mf2_opf_omp4kernel.cpp"
#include "poisson_mf2_opbf_omp4kernel.cpp"
#include "poisson_mf2_bc_omp4kernel.cpp"
#include "poisson_mf_edges_omp4kernel.cpp"
#include "poisson_mf_bedges_omp4kernel.cpp"
#include "poisson_mf_zero_omp4kernel.cpp"
#include "poisson_mf_bc_omp4kernel.cpp"
#include "poisson_mf2_apply_bc_vis_omp4kernel.cpp"
#include "poisson_mf2_apply_bc_omp4kernel.cpp"
#include "poisson_mf2_mass_omp4kernel.cpp"
#include "poisson_mf2_vis_omp4kernel.cpp"
#include "poisson_mf2_faces_vis_omp4kernel.cpp"
#include "poisson_mf2_omp4kernel.cpp"
#include "poisson_mf2_faces_omp4kernel.cpp"
#include "poisson_test_init_omp4kernel.cpp"
#include "poisson_test_bc_omp4kernel.cpp"
#include "poisson_test_error_omp4kernel.cpp"
#include "set_ic_omp4kernel.cpp"
#include "calc_dt_omp4kernel.cpp"
#include "advection_flux_omp4kernel.cpp"
#include "advection_faces_omp4kernel.cpp"
#include "advection_bc_omp4kernel.cpp"
#include "advection_numerical_flux_omp4kernel.cpp"
#include "advection_intermediate_vel_omp4kernel.cpp"
#include "pressure_bc_omp4kernel.cpp"
#include "pressure_bc2_omp4kernel.cpp"
#include "pressure_rhs_omp4kernel.cpp"
#include "pressure_update_vel_omp4kernel.cpp"
#include "viscosity_bc_omp4kernel.cpp"
#include "viscosity_rhs_omp4kernel.cpp"
#include "viscosity_reset_bc_omp4kernel.cpp"
#include "lift_drag_omp4kernel.cpp"
#include "save_values_omp4kernel.cpp"
#include "calc_h_omp4kernel.cpp"
#include "init_surface_omp4kernel.cpp"
#include "set_rkQ_omp4kernel.cpp"
#include "update_Q_omp4kernel.cpp"
#include "ls_advec_edges_omp4kernel.cpp"
#include "ls_advec_bedges_omp4kernel.cpp"
#include "ls_advec_flux_omp4kernel.cpp"
#include "ls_advec_rhs_omp4kernel.cpp"
#include "ls_sign_omp4kernel.cpp"
#include "ls_flux_omp4kernel.cpp"
#include "ls_bflux_omp4kernel.cpp"
#include "ls_copy_omp4kernel.cpp"
#include "ls_rhs_omp4kernel.cpp"
#include "ls_add_diff_omp4kernel.cpp"
#include "sigma_flux_omp4kernel.cpp"
#include "sigma_bflux_omp4kernel.cpp"
#include "sigma_mult_omp4kernel.cpp"
#include "diff_flux_omp4kernel.cpp"
#include "diff_bflux_omp4kernel.cpp"
#include "ls_reinit_check_omp4kernel.cpp"
#include "ls_step_omp4kernel.cpp"
