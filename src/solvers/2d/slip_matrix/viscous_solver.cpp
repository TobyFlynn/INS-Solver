#include "slip_matrix/2d/viscous_solver.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "dg_abort.h"

#include "timing.h"
#include "config.h"

extern Timing *timer;
extern Config *config;
extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

ViscousSolver::ViscousSolver(DGMesh2D *m) {
  u_bc = nullptr;
  v_bc = nullptr;
  nullspace = false;
  mesh = m;
  preconditioner = Preconditioners::NONE;
  max_iter = 5000;
  abs_tol = 1e-9;
  rel_tol = 1e-8;
  inv_mass_factor = 1.0;
  zero_input = false;
}

void ViscousSolver::set_matrix(ViscousMatrix2D *mat) {
  matrix = mat;
}

void ViscousSolver::set_bcs(op_dat u_bcs, op_dat v_bcs) {
  u_bc = u_bcs;
  v_bc = v_bcs;
}

void ViscousSolver::set_nullspace(bool ns) {
  nullspace = ns;
  if(nullspace) {
    dg_abort("Nullspace is not supported with ViscousSolver");
  }
}

void ViscousSolver::set_inv_mass_factor(DG_FP f) {
  inv_mass_factor = f;
}

void ViscousSolver::set_preconditioner(Preconditioners p) {
  preconditioner = p;
}

bool ViscousSolver::solve(op_dat u_rhs, op_dat v_rhs, op_dat u_ans, op_dat v_ans) {
  if(u_bc)
    matrix->apply_bc(u_rhs, v_rhs, u_bc, v_bc);

  if(zero_input) {
    op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                op_arg_dat(u_ans, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                op_arg_dat(v_ans, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }

  if(preconditioner == Preconditioners::NONE) {
   return conjugate_gradient(u_rhs, v_rhs, u_ans, v_ans);
  } else {
    return preconditioned_conjugate_gradient(u_rhs, v_rhs, u_ans, v_ans);
  }
}

void ViscousSolver::set_tol_and_iter(const double rtol, const double atol, const int maxiter) {
  max_iter = maxiter;
  abs_tol = atol;
  rel_tol = rtol;
}

bool ViscousSolver::conjugate_gradient(op_dat u_rhs, op_dat v_rhs, op_dat u_ans, op_dat v_ans) {
  DGTempDat cg_r[2];
  cg_r[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  cg_r[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat cg_p[2];
  cg_p[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  cg_p[1] = dg_dat_pool->requestTempDatCells(DG_NP);

  DG_FP residual = calc_res_explicit(u_rhs, v_rhs, u_ans, v_ans, cg_r[0].dat, cg_r[1].dat);

  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cg_p[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cg_p[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  
  DGTempDat cg_tmp[2];
  cg_tmp[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  cg_tmp[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  int iter = 0;
  // TODO add in relative tolerance
  while(residual > abs_tol && iter < max_iter) {
    // Calculate alpha
    matrix->mult(cg_p[0].dat, cg_p[1].dat, cg_tmp[0].dat, cg_tmp[1].dat);
    DG_FP tmp_alpha_0 = 0.0;
    DG_FP tmp_alpha_1 = 0.0;
    op_par_loop(cg_alpha, "cg_alpha", mesh->cells,
                op_arg_gbl(&tmp_alpha_0, 1, DG_FP_STR, OP_INC),
                op_arg_gbl(&tmp_alpha_1, 1, DG_FP_STR, OP_INC),
                op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_p[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_p[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
    DG_FP alpha = tmp_alpha_0 / tmp_alpha_1;

    // Update ans
    update_ans(alpha, cg_p[0].dat, cg_p[1].dat, u_ans, v_ans);
    // Compute denominator of beta now to avoid more temporary op_dats
    DG_FP tmp_beta_1 = 0.0;
    op_par_loop(cg_compute_beta, "cg_compute_beta", mesh->cells,
                op_arg_gbl(&tmp_beta_1, 1, DG_FP_STR, OP_INC),
                op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
    // Update r
    // residual = update_res(alpha, cg_tmp[0].dat, cg_tmp[1].dat, cg_r[0].dat, cg_r[1].dat);
    residual = calc_res_explicit(u_rhs, v_rhs, u_ans, v_ans, cg_r[0].dat, cg_r[1].dat);
    // Compute beta
    DG_FP tmp_beta_0 = 0.0;
    op_par_loop(cg_compute_beta, "cg_compute_beta", mesh->cells,
                op_arg_gbl(&tmp_beta_0, 1, DG_FP_STR, OP_INC),
                op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
    DG_FP beta = tmp_beta_0 / tmp_beta_1;
    // Update p
    update_p(beta, cg_r[0].dat, cg_r[1].dat, cg_p[0].dat, cg_p[1].dat);

    iter++;
  }

  op_printf("Final residual: %g\t\tIterations: %d\n", residual, iter);

  dg_dat_pool->releaseTempDatCells(cg_r[0]);
  dg_dat_pool->releaseTempDatCells(cg_r[1]);
  dg_dat_pool->releaseTempDatCells(cg_p[0]);
  dg_dat_pool->releaseTempDatCells(cg_p[1]);
  dg_dat_pool->releaseTempDatCells(cg_tmp[0]);
  dg_dat_pool->releaseTempDatCells(cg_tmp[1]);
  return residual < abs_tol;
}

// Actuall flexible CG now
bool ViscousSolver::preconditioned_conjugate_gradient(op_dat u_rhs, op_dat v_rhs, op_dat u_ans, op_dat v_ans) {
  DGTempDat cg_r[2];
  cg_r[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  cg_r[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat cg_z[2];
  cg_z[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  cg_z[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat cg_z_1[2];
  cg_z_1[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  cg_z_1[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat cg_p[2];
  cg_p[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  cg_p[1] = dg_dat_pool->requestTempDatCells(DG_NP);

  DG_FP residual = calc_res_explicit(u_rhs, v_rhs, u_ans, v_ans, cg_r[0].dat, cg_r[1].dat);

  precondition(cg_r[0].dat, cg_r[1].dat, cg_z[0].dat, cg_z[1].dat);

  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cg_p[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(cg_p[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  
  DGTempDat cg_tmp[2];
  cg_tmp[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  cg_tmp[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  int iter = 0;
  // TODO add in relative tolerance
  while(residual > abs_tol && iter < max_iter) {
    // Calculate alpha
    matrix->mult(cg_p[0].dat, cg_p[1].dat, cg_tmp[0].dat, cg_tmp[1].dat);
    DG_FP tmp_alpha_0 = 0.0;
    DG_FP tmp_alpha_1 = 0.0;
    op_par_loop(cg_alpha_pre, "cg_alpha_pre", mesh->cells,
                op_arg_gbl(&tmp_alpha_0, 1, DG_FP_STR, OP_INC),
                op_arg_gbl(&tmp_alpha_1, 1, DG_FP_STR, OP_INC),
                op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_p[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_p[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
    DG_FP alpha = tmp_alpha_0 / tmp_alpha_1;
    // Update ans
    update_ans(alpha, cg_p[0].dat, cg_p[1].dat, u_ans, v_ans);
    // Compute denominator of beta now to avoid more temporary op_dats
    DG_FP tmp_beta_1 = 0.0;
    op_par_loop(cg_compute_beta_pre, "cg_compute_beta_pre", mesh->cells,
                op_arg_gbl(&tmp_beta_1, 1, DG_FP_STR, OP_INC),
                op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
    // Update r
    residual = update_res(alpha, cg_tmp[0].dat, cg_tmp[1].dat, cg_r[0].dat, cg_r[1].dat);
    // residual = calc_res_explicit(u_rhs, v_rhs, u_ans, v_ans, cg_r[0].dat, cg_r[1].dat);
    // Save old z
    op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
                op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z_1[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
                op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z_1[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    // Update z
    precondition(cg_r[0].dat, cg_r[1].dat, cg_z[0].dat, cg_z[1].dat);
    // Compute beta
    DG_FP tmp_beta_0 = 0.0;   
    op_par_loop(cg_compute_beta_pre_flex, "cg_compute_beta_pre_flex", mesh->cells,
                op_arg_gbl(&tmp_beta_0, 1, DG_FP_STR, OP_INC),
                op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z_1[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z_1[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
/*
    op_par_loop(cg_compute_beta_pre, "cg_compute_beta_pre", mesh->cells,
                op_arg_gbl(&tmp_beta_0, 1, DG_FP_STR, OP_INC),
                op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
*/
    DG_FP beta = tmp_beta_0 / tmp_beta_1;
    // Update p
    update_p(beta, cg_z[0].dat, cg_z[1].dat, cg_p[0].dat, cg_p[1].dat);

    iter++;
  }

  op_printf("Final residual: %g\t\tIterations: %d\n", residual, iter);

  dg_dat_pool->releaseTempDatCells(cg_r[0]);
  dg_dat_pool->releaseTempDatCells(cg_r[1]);
  dg_dat_pool->releaseTempDatCells(cg_z[0]);
  dg_dat_pool->releaseTempDatCells(cg_z[1]);
  dg_dat_pool->releaseTempDatCells(cg_z_1[0]);
  dg_dat_pool->releaseTempDatCells(cg_z_1[1]);
  dg_dat_pool->releaseTempDatCells(cg_p[0]);
  dg_dat_pool->releaseTempDatCells(cg_p[1]);
  dg_dat_pool->releaseTempDatCells(cg_tmp[0]);
  dg_dat_pool->releaseTempDatCells(cg_tmp[1]);
  return residual < abs_tol;
}

// Helper functions
void ViscousSolver::update_ans(DG_FP alpha, op_dat p0, op_dat p1, op_dat u, op_dat v) {
  op_par_loop(cg_update_ans, "cg_update_ans", mesh->cells,
                op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
                op_arg_dat(p0, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(p1, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
}

DG_FP ViscousSolver::update_res(DG_FP alpha, op_dat Ap0, op_dat Ap1, op_dat r0, op_dat r1) {
  DG_FP residual = 0.0;
  op_par_loop(cg_update_res, "cg_update_res", mesh->cells,
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&residual, 1, DG_FP_STR, OP_INC),
              op_arg_dat(Ap0, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(Ap1, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(r0, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(r1, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  return sqrt(residual);
}

DG_FP ViscousSolver::calc_res_explicit(op_dat u_rhs, op_dat v_rhs, op_dat u_ans, op_dat v_ans, op_dat r0, op_dat r1) {
  matrix->mult(u_ans, v_ans, r0, r1);

  DG_FP residual = 0.0;
  op_par_loop(cg_residual, "cg_residual", mesh->cells,
                op_arg_gbl(&residual, 1, DG_FP_STR, OP_INC),
                op_arg_dat(u_rhs, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v_rhs, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(r0, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(r1, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  
  return sqrt(residual);
}

void ViscousSolver::update_p(DG_FP beta, op_dat r0, op_dat r1, op_dat p0, op_dat p1) {
  op_par_loop(cg_update_p, "cg_update_p", mesh->cells,
                op_arg_gbl(&beta, 1, DG_FP_STR, OP_READ),
                op_arg_dat(r0, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(r1, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(p0, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(p1, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
}

// Preconditioner functions
void ViscousSolver::precondition(op_dat u_res, op_dat v_res, op_dat u, op_dat v) {
  switch(preconditioner) {
    case Preconditioners::INV_MASS:
      pre_inv_mass(u_res, v_res, u, v);
      break;
    case Preconditioners::FACTOR_INV_MASS:
      pre_factor_inv_mass(u_res, v_res, u, v);
      break;
  }
}

void ViscousSolver::pre_inv_mass(op_dat u_res, op_dat v_res, op_dat u, op_dat v) {
  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(u_res, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  
  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(v_res, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  
  mesh->inv_mass(u);
  mesh->inv_mass(v);
}

void ViscousSolver::pre_factor_inv_mass(op_dat u_res, op_dat v_res, op_dat u, op_dat v) {
  pre_inv_mass(u_res, v_res, u, v);

  op_par_loop(constant_mult, "constant_mult", mesh->cells,
              op_arg_gbl(&inv_mass_factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  
  op_par_loop(constant_mult, "constant_mult", mesh->cells,
              op_arg_gbl(&inv_mass_factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}