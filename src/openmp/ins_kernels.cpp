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
extern double cubW[46];
extern double cubV[690];
extern double cubVDr[690];
extern double cubVDs[690];
extern double gF0Dr[105];
extern double gF0Ds[105];
extern double gF1Dr[105];
extern double gF1Ds[105];
extern double gF2Dr[105];
extern double gF2Ds[105];
extern double gaussW[7];
extern double gFInterp0[105];
extern double gFInterp1[105];
extern double gFInterp2[105];
extern double gF0DrR[105];
extern double gF0DsR[105];
extern double gF1DrR[105];
extern double gF1DsR[105];
extern double gF2DrR[105];
extern double gF2DsR[105];
extern double gFInterp0R[105];
extern double gFInterp1R[105];
extern double gFInterp2R[105];
extern double lift_drag_vec[5];

// header
#include "op_lib_cpp.h"

// user kernel files
#include "set_ic_kernel.cpp"
#include "calc_dt_kernel.cpp"
#include "advection_flux_kernel.cpp"
#include "advection_faces_kernel.cpp"
#include "advection_bc_kernel.cpp"
#include "advection_numerical_flux_kernel.cpp"
#include "advection_intermediate_vel_kernel.cpp"
#include "pressure_bc_kernel.cpp"
#include "pressure_rhs_kernel.cpp"
#include "pressure_update_vel_kernel.cpp"
#include "viscosity_bc_kernel.cpp"
#include "viscosity_rhs_kernel.cpp"
#include "viscosity_reset_bc_kernel.cpp"
#include "lift_drag_kernel.cpp"
#include "min_max_kernel.cpp"
#include "init_grid_kernel.cpp"
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
#include "cub_div_kernel.cpp"
#include "poisson_rhs_faces_kernel.cpp"
#include "poisson_rhs_bc_kernel.cpp"
#include "poisson_rhs_flux_kernel.cpp"
#include "poisson_rhs_J_kernel.cpp"
#include "poisson_rhs_qbc_kernel.cpp"
#include "poisson_rhs_qflux_kernel.cpp"
#include "poisson_test_init_kernel.cpp"
#include "poisson_test_bc_kernel.cpp"
#include "poisson_test_error_kernel.cpp"
