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

// header
#include "op_lib_cpp.h"

// user kernel files
#include "init_grid_seqkernel.cpp"
#include "set_ic_seqkernel.cpp"
#include "calc_dt_seqkernel.cpp"
#include "advection_flux_seqkernel.cpp"
#include "advection_faces_seqkernel.cpp"
#include "advection_bc_seqkernel.cpp"
#include "advection_numerical_flux_seqkernel.cpp"
#include "advection_intermediate_vel_seqkernel.cpp"
#include "pressure_bc_seqkernel.cpp"
#include "pressure_rhs_seqkernel.cpp"
#include "pressure_nBC_seqkernel.cpp"
#include "pressure_bc2_seqkernel.cpp"
#include "pressure_bc3_seqkernel.cpp"
#include "pressure_update_vel_seqkernel.cpp"
#include "reset_nBC_seqkernel.cpp"
#include "viscosity_faces_seqkernel.cpp"
#include "viscosity_rhs_seqkernel.cpp"
#include "viscosity_bc_seqkernel.cpp"
#include "viscosity_reset_bc_seqkernel.cpp"
#include "setup_poisson_seqkernel.cpp"
#include "set_tau_seqkernel.cpp"
#include "set_tau_bc_seqkernel.cpp"
#include "poisson_rhs_faces_seqkernel.cpp"
#include "poisson_rhs_bc_seqkernel.cpp"
#include "poisson_rhs_du_seqkernel.cpp"
#include "poisson_rhs_qbc_seqkernel.cpp"
#include "poisson_rhs_fluxq_seqkernel.cpp"
#include "poisson_rhs_J_seqkernel.cpp"
#include "div_seqkernel.cpp"
#include "curl_seqkernel.cpp"
#include "grad_seqkernel.cpp"
#include "poisson_test_init_seqkernel.cpp"
#include "poisson_test_bc_seqkernel.cpp"
#include "poisson_test_set_rhs_seqkernel.cpp"
#include "poisson_test_error_seqkernel.cpp"
