//
// auto-generated by op2.py
//

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
#include "init_nodes_seqkernel.cpp"
#include "init_grid_seqkernel.cpp"
#include "init_cubature_grad_seqkernel.cpp"
#include "init_cubature_seqkernel.cpp"
#include "init_cubature_OP_seqkernel.cpp"
#include "gauss_reverse_seqkernel.cpp"
#include "init_gauss_seqkernel.cpp"
#include "gauss_tau_seqkernel.cpp"
#include "gauss_tau_bc_seqkernel.cpp"
#include "init_gauss_grad_seqkernel.cpp"
#include "init_gauss_grad2_seqkernel.cpp"
#include "init_gauss_grad_neighbour_seqkernel.cpp"
#include "gauss_grad_faces_seqkernel.cpp"
#include "gauss_op_seqkernel.cpp"
#include "gauss_gfi_faces_seqkernel.cpp"
#include "div_seqkernel.cpp"
#include "curl_seqkernel.cpp"
#include "grad_seqkernel.cpp"
#include "glb_ind_kernel_seqkernel.cpp"
#include "glb_ind_kernelBC_seqkernel.cpp"
#include "poisson_mf2_op_seqkernel.cpp"
#include "poisson_mf2_opf_seqkernel.cpp"
#include "poisson_mf2_opbf_seqkernel.cpp"
#include "poisson_mf2_bc_seqkernel.cpp"
#include "poisson_mf2_apply_bc_seqkernel.cpp"
#include "poisson_mf2_mass_seqkernel.cpp"
#include "poisson_mf2_seqkernel.cpp"
#include "poisson_mf2_faces_seqkernel.cpp"
#include "poisson_test_init_seqkernel.cpp"
#include "poisson_test_bc_seqkernel.cpp"
#include "poisson_test_error_seqkernel.cpp"
#include "set_ic_seqkernel.cpp"
#include "calc_dt_seqkernel.cpp"
#include "advection_flux_seqkernel.cpp"
#include "advection_faces_seqkernel.cpp"
#include "advection_bc_seqkernel.cpp"
#include "advection_numerical_flux_seqkernel.cpp"
#include "advection_intermediate_vel_seqkernel.cpp"
#include "pressure_bc_seqkernel.cpp"
#include "pressure_bc2_seqkernel.cpp"
#include "pressure_rhs_seqkernel.cpp"
#include "pressure_update_vel_seqkernel.cpp"
#include "viscosity_bc_seqkernel.cpp"
#include "viscosity_rhs_seqkernel.cpp"
#include "viscosity_reset_bc_seqkernel.cpp"
#include "lift_drag_seqkernel.cpp"
#include "save_values_seqkernel.cpp"
#include "init_surface_seqkernel.cpp"
#include "set_rkQ_seqkernel.cpp"
#include "update_Q_seqkernel.cpp"
#include "ls_advec_edges_seqkernel.cpp"
#include "ls_advec_bedges_seqkernel.cpp"
#include "ls_advec_flux_seqkernel.cpp"
#include "ls_advec_rhs_seqkernel.cpp"
