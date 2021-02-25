//
// auto-generated by op2.py
//

// global constants
double gam_ompkernel;
double mu_ompkernel;
double nu_ompkernel;
double bc_mach_ompkernel;
double bc_alpha_ompkernel;
double bc_p_ompkernel;
double bc_u_ompkernel;
double bc_v_ompkernel;
int FMASK_ompkernel[15];
double ic_u_ompkernel;
double ic_v_ompkernel;

// header
#include "op_lib_cpp.h"

void op_decl_const_char(int dim, char const *type,
  int size, char *dat, char const *name){
  if(!strcmp(name, "gam")) {
    memcpy(&gam_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gam_ompkernel)
  } else if(!strcmp(name, "mu")) {
    memcpy(&mu_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:mu_ompkernel)
  } else if(!strcmp(name, "nu")) {
    memcpy(&nu_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:nu_ompkernel)
  } else if(!strcmp(name, "bc_mach")) {
    memcpy(&bc_mach_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_mach_ompkernel)
  } else if(!strcmp(name, "bc_alpha")) {
    memcpy(&bc_alpha_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_alpha_ompkernel)
  } else if(!strcmp(name, "bc_p")) {
    memcpy(&bc_p_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_p_ompkernel)
  } else if(!strcmp(name, "bc_u")) {
    memcpy(&bc_u_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_u_ompkernel)
  } else if(!strcmp(name, "bc_v")) {
    memcpy(&bc_v_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_v_ompkernel)
  } else if(!strcmp(name, "FMASK")) {
    memcpy(FMASK_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:FMASK_ompkernel[:15])
  } else if(!strcmp(name, "ic_u")) {
    memcpy(&ic_u_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:ic_u_ompkernel)
  } else if(!strcmp(name, "ic_v")) {
    memcpy(&ic_v_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:ic_v_ompkernel)
  }
}
// user kernel files
#include "init_grid_omp4kernel_func.cpp"
#include "set_ic_omp4kernel_func.cpp"
#include "calc_dt_omp4kernel_func.cpp"
#include "advection_flux_omp4kernel_func.cpp"
#include "advection_faces_omp4kernel_func.cpp"
#include "advection_bc_omp4kernel_func.cpp"
#include "advection_numerical_flux_omp4kernel_func.cpp"
#include "advection_intermediate_vel_omp4kernel_func.cpp"
#include "pressure_bc_omp4kernel_func.cpp"
#include "pressure_rhs_omp4kernel_func.cpp"
#include "pressure_bc2_omp4kernel_func.cpp"
#include "pressure_bc3_omp4kernel_func.cpp"
#include "pressure_update_vel_omp4kernel_func.cpp"
#include "viscosity_faces_omp4kernel_func.cpp"
#include "viscosity_rhs_omp4kernel_func.cpp"
#include "viscosity_bc_omp4kernel_func.cpp"
#include "viscosity_reset_bc_omp4kernel_func.cpp"
#include "setup_poisson_omp4kernel_func.cpp"
#include "set_tau_omp4kernel_func.cpp"
#include "set_tau_bc_omp4kernel_func.cpp"
#include "poisson_rhs_faces_omp4kernel_func.cpp"
#include "poisson_rhs_bc_omp4kernel_func.cpp"
#include "poisson_rhs_du_omp4kernel_func.cpp"
#include "poisson_rhs_qbc_omp4kernel_func.cpp"
#include "poisson_rhs_fluxq_omp4kernel_func.cpp"
#include "poisson_rhs_J_omp4kernel_func.cpp"
#include "div_omp4kernel_func.cpp"
#include "curl_omp4kernel_func.cpp"
#include "grad_omp4kernel_func.cpp"
#include "poisson_test_init_omp4kernel_func.cpp"
#include "poisson_test_bc_omp4kernel_func.cpp"
#include "poisson_test_set_rhs_omp4kernel_func.cpp"
#include "poisson_test_error_omp4kernel_func.cpp"
