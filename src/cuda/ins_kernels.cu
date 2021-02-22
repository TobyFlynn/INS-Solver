//
// auto-generated by op2.py
//

//global constants
#ifndef MAX_CONST_SIZE
#define MAX_CONST_SIZE 128
#endif

__constant__ double gam_cuda;
__constant__ double mu_cuda;
__constant__ double nu_cuda;
__constant__ double bc_mach_cuda;
__constant__ double bc_alpha_cuda;
__constant__ double bc_p_cuda;
__constant__ double bc_u_cuda;
__constant__ double bc_v_cuda;
__constant__ int FMASK_cuda[15];
__constant__ double ic_u_cuda;
__constant__ double ic_v_cuda;

//header
#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

void op_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){
  if (!OP_hybrid_gpu) return;
  if (!strcmp(name,"gam")) {
    cutilSafeCall(cudaMemcpyToSymbol(gam_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"mu")) {
    cutilSafeCall(cudaMemcpyToSymbol(mu_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"nu")) {
    cutilSafeCall(cudaMemcpyToSymbol(nu_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_mach")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_mach_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_alpha")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_alpha_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_p")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_p_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_u")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_u_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"bc_v")) {
    cutilSafeCall(cudaMemcpyToSymbol(bc_v_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"FMASK")) {
    cutilSafeCall(cudaMemcpyToSymbol(FMASK_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"ic_u")) {
    cutilSafeCall(cudaMemcpyToSymbol(ic_u_cuda, dat, dim*size));
  }
  else
  if (!strcmp(name,"ic_v")) {
    cutilSafeCall(cudaMemcpyToSymbol(ic_v_cuda, dat, dim*size));
  }
  else
  {
    printf("error: unknown const name\n"); exit(1);
  }
}

//user kernel files
#include "init_grid_kernel.cu"
#include "set_ic_kernel.cu"
#include "calc_dt_kernel.cu"
#include "advection_flux_kernel.cu"
#include "advection_faces_kernel.cu"
#include "advection_bc_kernel.cu"
#include "advection_numerical_flux_kernel.cu"
#include "advection_intermediate_vel_kernel.cu"
#include "pressure_bc_kernel.cu"
#include "pressure_rhs_kernel.cu"
#include "pressure_bc2_kernel.cu"
#include "pressure_bc3_kernel.cu"
#include "pressure_update_vel_kernel.cu"
#include "viscosity_faces_kernel.cu"
#include "viscosity_set_bc_kernel.cu"
#include "viscosity_rhs_kernel.cu"
#include "viscosity_bc_kernel.cu"
#include "setup_poisson_kernel.cu"
#include "set_tau_kernel.cu"
#include "set_tau_bc_kernel.cu"
#include "poisson_rhs_faces_kernel.cu"
#include "poisson_rhs_bc_kernel.cu"
#include "poisson_rhs_du_kernel.cu"
#include "poisson_rhs_qbc_kernel.cu"
#include "poisson_rhs_fluxq_kernel.cu"
#include "poisson_rhs_J_kernel.cu"
#include "div_kernel.cu"
#include "curl_kernel.cu"
#include "grad_kernel.cu"
#include "poisson_test_init_kernel.cu"
#include "poisson_test_bc_kernel.cu"
#include "poisson_test_set_rhs_kernel.cu"
#include "poisson_test_error_kernel.cu"
