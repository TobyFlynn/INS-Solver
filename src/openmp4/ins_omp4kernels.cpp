//
// auto-generated by op2.py
//

// header
#include "op_lib_cpp.h"

// user kernel files
#include "set_ic_omp4kernel.cpp"
#include "calc_dt_omp4kernel.cpp"
#include "advection_flux_omp4kernel.cpp"
#include "advection_faces_omp4kernel.cpp"
#include "advection_bc_omp4kernel.cpp"
#include "advection_numerical_flux_omp4kernel.cpp"
#include "advection_intermediate_vel_omp4kernel.cpp"
#include "pressure_bc_omp4kernel.cpp"
#include "pressure_rhs_omp4kernel.cpp"
#include "pressure_update_vel_omp4kernel.cpp"
#include "viscosity_bc_omp4kernel.cpp"
#include "viscosity_rhs_omp4kernel.cpp"
#include "viscosity_reset_bc_omp4kernel.cpp"
#include "lift_drag_omp4kernel.cpp"
#include "min_max_omp4kernel.cpp"
#include "init_grid_omp4kernel.cpp"
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
#include "poisson_test_init_omp4kernel.cpp"
#include "poisson_test_bc_omp4kernel.cpp"
#include "poisson_test_set_rhs_omp4kernel.cpp"
#include "poisson_test_error_omp4kernel.cpp"
