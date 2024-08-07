#include "solvers/2d/ins_solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>
#include <sstream>

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "dg_utils.h"
#include "dg_abort.h"

#include "timing.h"
#include "config.h"
#include "dg_linear_solvers/petsc_amg.h"
#include "dg_linear_solvers/petsc_block_jacobi.h"
#include "dg_linear_solvers/petsc_pmultigrid.h"
#include "dg_linear_solvers/initial_guess_extrapolation.h"

extern Timing *timer;
extern Config *config;
extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

INSSolverBase2D::INSSolverBase2D(DGMesh2D *m) {
  mesh = m;

  read_options();
  init_dats();

  std::string name;
  for(int i = 0; i < 2; i++) {
    name = "ins_solver_vel0" + std::to_string(i);
    vel[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "ins_solver_vel1" + std::to_string(i);
    vel[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "ins_solver_n0" + std::to_string(i);
    n[0][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "ins_solver_n1" + std::to_string(i);
    n[1][i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
  }
  pr = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_pr");

  dPdN[0] = op_decl_dat(mesh->cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_dPdN0");
  dPdN[1] = op_decl_dat(mesh->cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_dPdN1");

  prev_time = 0.0;
  time = 0.0;
  currentInd = 0;
}

INSSolverBase2D::INSSolverBase2D(DGMesh2D *m, const std::string &filename) {
  mesh = m;

  read_options();
  init_dats();

  vel[0][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel00");
  vel[1][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel10");
  vel[0][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel01");
  vel[1][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel11");
  n[0][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n00");
  n[1][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n10");
  n[0][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n01");
  n[1][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n11");
  pr = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_pr");
  dPdN[0] = op_decl_dat_hdf5(mesh->cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, filename.c_str(), "ins_solver_dPdN0");
  dPdN[1] = op_decl_dat_hdf5(mesh->cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, filename.c_str(), "ins_solver_dPdN1");

  op_get_const_hdf5("time", 1, "double", (char *)&time, filename.c_str());
  op_get_const_hdf5("prev_time", 1, "double", (char *)&prev_time, filename.c_str());
  op_get_const_hdf5("currentInd", 1, "int", (char *)&currentInd, filename.c_str());
  op_get_const_hdf5("it_pre_sub_cycle", 1, "int", (char *)&it_pre_sub_cycle, filename.c_str());
}

void INSSolverBase2D::read_options() {
  // Get vector of values that should be saved in output HDF5 file
  std::string save_tmp = "all";
  config->getStr("io", "values_to_save", save_tmp);
  if(save_tmp == "all")
    save_tmp = "velocity,pressure,non_linear,surface_tension,intermediate_velocities,mu,rho,level_set,curvature,art_vis";
  std::stringstream tmp_ss(save_tmp);
  std::string val_str;
  while(std::getline(tmp_ss, val_str, ',')) {
    values_to_save.insert(val_str);
  }

  int tmp_div = 0;
  config->getInt("solver-options", "div_div", tmp_div);
  pr_projection_method = tmp_div;
  config->getInt("solver-options", "sub_cycle", sub_cycles);
  it_pre_sub_cycle = 1;
  config->getInt("solver-options", "num_iter_before_sub_cycle", it_pre_sub_cycle);
  it_pre_sub_cycle = it_pre_sub_cycle > 1 ? it_pre_sub_cycle : 1;
  int tmp_eig = 1;
  config->getInt("solver-options", "extrapolate_initial_guess", tmp_eig);
  extrapolate_initial_guess = tmp_eig == 1;
  int tmp_shock = 0;
  config->getInt("solver-options", "shock_capturing", tmp_shock);
  shock_capturing = tmp_shock == 1;
  int tmp_oia = 0;
  config->getInt("solver-options", "over_int_advec", tmp_oia);
  over_int_advec = tmp_oia == 1;
  int tmp_grav = 0;
  config->getInt("solver-options", "gravity", tmp_grav);
  gravity = tmp_grav == 1;

  // if(gravity && sub_cycles > 0)
  //   dg_abort("Gravity not supported with subcycling currently");

  filter_alpha = 18.0;
  config->getDouble("filter", "alpha", filter_alpha);
  filter_sp = 32;
  config->getInt("filter", "sp", filter_sp);
  filter_Nc = 0;
  config->getInt("filter", "Nc", filter_Nc);

  shock_cap_max_diff = 0.0;
  config->getDouble("shock-capturing", "max_art_diff", shock_cap_max_diff);
  shock_cap_smooth_tol = 0.5;
  config->getDouble("shock-capturing", "smooth_tol", shock_cap_smooth_tol);
  shock_cap_discon_tol = 1.5;
  config->getDouble("shock-capturing", "discont_tol", shock_cap_discon_tol);


}

void INSSolverBase2D::init_dats() {
  std::string name;
  for(int i = 0; i < 2; i++) {
    name     = "ins_solver_velT" + std::to_string(i);
    velT[i]  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name     = "ins_solver_velTT" + std::to_string(i);
    velTT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
  }

  bc_types = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_bc_types");

  if(pr_projection_method == 1 || pr_projection_method == 2) {
    proj_h = op_decl_dat(mesh->cells, 1, DG_FP_STR, (DG_FP *)NULL, "ins_solver_proj_h");
  }

  bc_data = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_bc_data");

  if(shock_capturing) {
    nodes_data        = op_decl_dat(mesh->nodes, 1, DG_FP_STR, (DG_FP *)NULL, "ins_solver_nodes_data");
    nodes_count       = op_decl_dat(mesh->nodes, 1, DG_FP_STR, (DG_FP *)NULL, "ins_solver_nodes_count");
  }
}

void INSSolverBase2D::setup_pressure_viscous_solvers(LinearSolver *pr_solver, LinearSolver *vis_solver) {
  int pr_tmp = 0;
  config->getInt("pressure-solve", "nullspace", pr_tmp);
  pr_solver->set_nullspace(pr_tmp == 1);
  double r_tol, a_tol;
  if(std::is_same<DG_FP,double>::value) {
    r_tol = 1e-8;
    a_tol = 1e-9;
  } else {
    r_tol = 1e-5;
    a_tol = 1e-6;
  }
  int max_iter = 500;
  config->getDouble("pressure-solve", "r_tol", r_tol);
  config->getDouble("pressure-solve", "a_tol", a_tol);
  config->getInt("pressure-solve", "max_iter", max_iter);
  pr_solver->set_tol_and_iter(r_tol, a_tol, max_iter);

  if(vis_solver) {
    int vis_tmp = 0;
    config->getInt("viscous-solve", "nullspace", vis_tmp);
    vis_solver->set_nullspace(vis_tmp == 1);
    if(std::is_same<DG_FP,double>::value) {
      r_tol = 1e-8;
      a_tol = 1e-9;
    } else {
      r_tol = 1e-5;
      a_tol = 1e-6;
    }
    max_iter = 5000;
    config->getDouble("viscous-solve", "r_tol", r_tol);
    config->getDouble("viscous-solve", "a_tol", a_tol);
    config->getInt("viscous-solve", "max_iter", max_iter);
    vis_solver->set_tol_and_iter(r_tol, a_tol, max_iter);
  }
}

INSSolverBase2D::~INSSolverBase2D() {

}

void INSSolverBase2D::init() {
  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  shock_cap_max_diff *= h;
  op_printf("h: %g\n", h);

  if(pr_projection_method == 1 || pr_projection_method == 2) {
    DGTempDat tmp_npf = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
    zero_dat(tmp_npf.dat);

    op_par_loop(ins_proj_setup_0, "ins_proj_setup_0", mesh->faces,
                op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
                op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));

    op_par_loop(ins_proj_setup_1, "ins_proj_setup_1", mesh->cells,
                op_arg_dat(tmp_npf.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));
    dg_dat_pool->releaseTempDatCells(tmp_npf);
  }

  // Calculate area of the domain
  DGTempDat tmp_dat = dg_dat_pool->requestTempDatCells(DG_NP);
  DG_FP value = 1.0;
  op_par_loop(set_val_dg_np, "set_val_dg_np", mesh->cells,
              op_arg_gbl(&value, 1, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_dat.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  mesh->mass(tmp_dat.dat);
  domain_area = 0.0;
  op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
              op_arg_gbl(&domain_area, 1, DG_FP_STR, OP_INC),
              op_arg_dat(tmp_dat.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  dg_dat_pool->releaseTempDatCells(tmp_dat);
}


void INSSolverBase2D::advec_current_non_linear() {
  DGTempDat f[2][2];
  f[0][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[0][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  // Calculate flux values
  op_par_loop(ins_advec_flux_2d, "ins_advec_flux_2d", mesh->cells,
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div(f[0][0].dat, f[0][1].dat, n[currentInd][0]);
  mesh->div(f[1][0].dat, f[1][1].dat, n[currentInd][1]);

  dg_dat_pool->releaseTempDatCells(f[0][0]);
  dg_dat_pool->releaseTempDatCells(f[0][1]);
  dg_dat_pool->releaseTempDatCells(f[1][0]);
  dg_dat_pool->releaseTempDatCells(f[1][1]);

  DGTempDat tmp_advec_flux0 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_advec_flux1 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  zero_dat(tmp_advec_flux0.dat);
  zero_dat(tmp_advec_flux1.dat);

  op_par_loop(ins_advec_faces_2d, "ins_advec_faces_2d", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale,  -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0],  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1],  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_advec_flux0.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux1.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  // Enforce BCs
  if(mesh->bface2cells) {
    op_par_loop(ins_advec_bc_2d, "ins_advec_bc_2d", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_advec_flux0.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_advec_flux1.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux0.dat, 1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux1.dat, 1.0, n[currentInd][1]);

  dg_dat_pool->releaseTempDatCells(tmp_advec_flux0);
  dg_dat_pool->releaseTempDatCells(tmp_advec_flux1);
}

void INSSolverBase2D::advec_current_non_linear_over_int() {
  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int");
  DGTempDat f[2][2];
  f[0][0] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  f[0][1] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  f[1][0] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  f[1][1] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);

  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int - Interp");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, vel[currentInd][0], 0.0, f[0][0].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, vel[currentInd][1], 0.0, f[0][1].dat);
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int - Interp");

  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int - 0");
  op_par_loop(ins_2d_advec_oi_0, "ins_2d_advec_oi_0", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(f[0][1].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(f[1][0].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int - 0");

  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int - div");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, f[0][0].dat, 0.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, f[0][1].dat, 1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, f[1][0].dat, 0.0, n[currentInd][1]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, f[1][1].dat, 1.0, n[currentInd][1]);
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int - div");

  dg_dat_pool->releaseTempDatCells(f[0][0]);
  dg_dat_pool->releaseTempDatCells(f[0][1]);
  dg_dat_pool->releaseTempDatCells(f[1][0]);
  dg_dat_pool->releaseTempDatCells(f[1][1]);

  DGTempDat fU_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  DGTempDat fV_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);

  zero_dat(fU_cub.dat);
  zero_dat(fV_cub.dat);

  op_par_loop(ins_2d_advec_oi_1, "ins_2d_advec_oi_1", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fU_cub.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(fV_cub.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE));

  if(mesh->bface2cells) {
    op_par_loop(ins_2d_advec_oi_2, "ins_2d_advec_oi_2", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(fU_cub.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC),
                op_arg_dat(fV_cub.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, fU_cub.dat, -1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, fV_cub.dat, -1.0, n[currentInd][1]);

  dg_dat_pool->releaseTempDatCells(fU_cub);
  dg_dat_pool->releaseTempDatCells(fV_cub);
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int");
}

void INSSolverBase2D::advec_standard() {
  if(over_int_advec) {
    advec_current_non_linear_over_int();
  } else {
    advec_current_non_linear();
  }

  // Calculate the intermediate velocity values
  if(gravity) {
    op_par_loop(ins_advec_intermediate_vel_grav_2d, "ins_advec_intermediate_vel_grav_2d", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  } else {
    op_par_loop(ins_advec_intermediate_vel_2d, "ins_advec_intermediate_vel_2d", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }
}

void INSSolverBase2D::advec_standard(op_dat fx, op_dat fy, op_dat fx_old, op_dat fy_old) {
  if(over_int_advec) {
    advec_current_non_linear_over_int();
  } else {
    advec_current_non_linear();
  }

  // Calculate the intermediate velocity values
  if(gravity) {
    op_par_loop(ins_advec_intermediate_vel_force_grav_2d, "ins_advec_intermediate_vel_force_grav_2d", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fx_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  } else {
    op_par_loop(ins_advec_intermediate_vel_force_2d, "ins_advec_intermediate_vel_force_2d", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fx_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }
}

void INSSolverBase2D::advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat u_out,
                                      op_dat v_out, const double t) {
  double t0 = prev_time;
  double t1 = time;
  double tI = t;
  DGTempDat advec_sc[2];
  advec_sc[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc[1] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_advec_sc_interp_2d, "ins_advec_sc_interp_2d", mesh->cells,
              op_arg_gbl(&t0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&t1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&tI, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  DGTempDat f[2][2];
  f[0][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[0][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][1] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_advec_sc_rhs_0_2d, "ins_advec_sc_rhs_0_2d", mesh->cells,
              op_arg_dat(u_in, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v_in, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div(f[0][0].dat, f[0][1].dat, u_out);
  mesh->div(f[1][0].dat, f[1][1].dat, v_out);

  dg_dat_pool->releaseTempDatCells(f[0][0]);
  dg_dat_pool->releaseTempDatCells(f[0][1]);
  dg_dat_pool->releaseTempDatCells(f[1][0]);
  dg_dat_pool->releaseTempDatCells(f[1][1]);

  DGTempDat tmp_advec_flux0 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_advec_flux1 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  zero_dat(tmp_advec_flux0.dat);
  zero_dat(tmp_advec_flux1.dat);

  op_par_loop(ins_advec_sc_rhs_1_2d, "ins_advec_sc_rhs_1_2d", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny,      -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale,  -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(u_in,  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v_in,  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[1].dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_advec_flux0.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux1.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  if(mesh->bface2cells) {
    op_par_loop(ins_advec_sc_rhs_2_2d, "ins_advec_sc_rhs_2_2d", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc[0].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc[1].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_advec_flux0.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_advec_flux1.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

  dg_dat_pool->releaseTempDatCells(advec_sc[0]);
  dg_dat_pool->releaseTempDatCells(advec_sc[1]);

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux0.dat, 1.0, u_out);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux1.dat, 1.0, v_out);

  dg_dat_pool->releaseTempDatCells(tmp_advec_flux0);
  dg_dat_pool->releaseTempDatCells(tmp_advec_flux1);
}

void INSSolverBase2D::advec_sub_cycle_rhs_over_int(op_dat u_in, op_dat v_in,
                                  op_dat u_out, op_dat v_out, const double t) {
  double t0 = prev_time;
  double t1 = time;
  double tI = t;
  DGTempDat advec_sc[2];
  advec_sc[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc[1] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_advec_sc_interp_2d, "ins_advec_sc_interp_2d", mesh->cells,
              op_arg_gbl(&t0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&t1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&tI, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  DGTempDat f[2][2];
  f[0][0] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  f[0][1] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  f[1][0] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  f[1][1] = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, u_in, 0.0, f[0][0].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, v_in, 0.0, f[0][1].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, advec_sc[0].dat, 0.0, f[1][0].dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, advec_sc[1].dat, 0.0, f[1][1].dat);

  op_par_loop(ins_2d_advec_sc_rhs_oi_0, "ins_2d_advec_sc_rhs_oi_0", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(f[0][1].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(f[1][0].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(f[1][1].dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, f[0][0].dat, 0.0, u_out);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, f[0][1].dat, 1.0, u_out);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, f[1][0].dat, 0.0, v_out);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, f[1][1].dat, 1.0, v_out);

  dg_dat_pool->releaseTempDatCells(f[0][0]);
  dg_dat_pool->releaseTempDatCells(f[0][1]);
  dg_dat_pool->releaseTempDatCells(f[1][0]);
  dg_dat_pool->releaseTempDatCells(f[1][1]);

  DGTempDat fU_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  DGTempDat fV_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);

  zero_dat(fU_cub.dat);
  zero_dat(fV_cub.dat);

  op_par_loop(ins_2d_advec_sc_rhs_oi_1, "ins_2d_advec_sc_rhs_oi_1", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(u_in, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v_in, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[1].dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fU_cub.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(fV_cub.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE));

  if(mesh->bface2cells) {
    op_par_loop(ins_2d_advec_sc_rhs_oi_2, "ins_2d_advec_sc_rhs_oi_2", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc[0].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc[1].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(fU_cub.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC),
                op_arg_dat(fV_cub.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, fU_cub.dat, -1.0, u_out);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, fV_cub.dat, -1.0, v_out);

  dg_dat_pool->releaseTempDatCells(fU_cub);
  dg_dat_pool->releaseTempDatCells(fV_cub);
  dg_dat_pool->releaseTempDatCells(advec_sc[0]);
  dg_dat_pool->releaseTempDatCells(advec_sc[1]);
}

void INSSolverBase2D::advec_sub_cycle_rk_step(const DG_FP time_sc, const DG_FP rk_dt, op_dat u, op_dat v) {
  // Request temporary dats
  DGTempDat advec_sc_rk[3][2];
  advec_sc_rk[0][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[0][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[1][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[1][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[2][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[2][1] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_advec_sc_copy_2d, "ins_advec_sc_copy_2d", mesh->cells,
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc_rk[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  for(int rk_step = 0; rk_step < 3; rk_step++) {
    double rk_time = time_sc;
    if(rk_step == 1) rk_time += rk_dt;
    if(rk_step == 2) rk_time += 0.5 * rk_dt;
    const int rk_ind = rk_step == 2 ? 2 : rk_step + 1;
    if(over_int_advec) {
      advec_sub_cycle_rhs_over_int(advec_sc_rk[0][0].dat, advec_sc_rk[0][1].dat,
                          advec_sc_rk[rk_ind][0].dat, advec_sc_rk[rk_ind][1].dat,
                          rk_time);
    } else {
      advec_sub_cycle_rhs(advec_sc_rk[0][0].dat, advec_sc_rk[0][1].dat,
                          advec_sc_rk[rk_ind][0].dat, advec_sc_rk[rk_ind][1].dat,
                          rk_time);
    }
    // Set up next step
    if(rk_step == 0) {
      op_par_loop(ins_advec_sc_rk_0_2d, "ins_advec_sc_rk_0_2d", mesh->cells,
                  op_arg_gbl(&rk_dt, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                  op_arg_dat(advec_sc_rk[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    } else if(rk_step == 1) {
      op_par_loop(ins_advec_sc_rk_1_2d, "ins_advec_sc_rk_1_2d", mesh->cells,
                  op_arg_gbl(&rk_dt, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(advec_sc_rk[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(advec_sc_rk[2][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[2][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                  op_arg_dat(advec_sc_rk[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    }
  }
  // Update velT
  op_par_loop(ins_advec_sc_rk_2_2d, "ins_advec_sc_rk_2_2d", mesh->cells,
              op_arg_gbl(&rk_dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[2][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[2][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  dg_dat_pool->releaseTempDatCells(advec_sc_rk[0][0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[0][1]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[1][0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[1][1]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[2][0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[2][1]);
}

void INSSolverBase2D::advec_sub_cycle() {
  op_par_loop(ins_advec_sc_copy_2d, "ins_advec_sc_copy_2d", mesh->cells,
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  // Advance 2 * number of subcycles
  const int num_steps = std::ceil((time + dt - prev_time) / sub_cycle_dt);
  const DG_FP tmp_sub_cycle_dt = (time + dt - prev_time) / (DG_FP)num_steps;
  for(int i = 0; i < num_steps; i++) {
    advec_sub_cycle_rk_step(prev_time + i * tmp_sub_cycle_dt, tmp_sub_cycle_dt, velT[0], velT[1]);
  }

  DGTempDat advec_sc_tmp[2];
  advec_sc_tmp[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_tmp[1] = dg_dat_pool->requestTempDatCells(DG_NP);

  // Other velocity
  op_par_loop(ins_advec_sc_copy_2d, "ins_advec_sc_copy_2d", mesh->cells,
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  // Advance number of subcycles
  for(int i = 0; i < sub_cycles; i++) {
    advec_sub_cycle_rk_step(time + i * sub_cycle_dt, sub_cycle_dt, advec_sc_tmp[0].dat, advec_sc_tmp[1].dat);
  }

  // Get final velT
  if(gravity) {
    op_par_loop(ins_advec_sc_update_grav_2d, "ins_advec_sc_update_grav_2d", mesh->cells,
                op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  } else {
    op_par_loop(ins_advec_sc_update_2d, "ins_advec_sc_update_2d", mesh->cells,
                op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  }

  dg_dat_pool->releaseTempDatCells(advec_sc_tmp[0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_tmp[1]);

  // Needed for Pressure boundary conditions
  if(over_int_advec) {
    advec_current_non_linear_over_int();
  } else {
    advec_current_non_linear();
  }
}

void INSSolverBase2D::advec_sub_cycle(op_dat fx, op_dat fy, op_dat fx_old, op_dat fy_old) {
  advec_sub_cycle();

  // Calculate the intermediate velocity values
  if(gravity) {
    op_par_loop(ins_2d_advec_sc_force_grav_2d, "ins_2d_advec_sc_force_grav_2d", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(fx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fx_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  } else {
    op_par_loop(ins_2d_advec_sc_force_2d, "ins_2d_advec_sc_force_2d", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(fx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fx_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fy_old, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  }
}

void INSSolverBase2D::project_velocity_mat_mult(op_dat u, op_dat v,
          op_dat u_out, op_dat v_out, op_dat pen, op_dat pen_f) {
  timer->startTimer("pr_proj - mult");
  DGTempDat div_tmp = dg_dat_pool->requestTempDatCells(DG_NP);

  timer->startTimer("pr_proj - mult - div");
  mesh->div(u, v, div_tmp.dat);
  timer->endTimer("pr_proj - mult - div");

  timer->startTimer("pr_proj - mult - mass");
  mesh->mass(div_tmp.dat);
  timer->endTimer("pr_proj - mult - mass");

  timer->startTimer("pr_proj - mult - drst");
  op2_gemv(mesh, true, 1.0, DGConstants::DR, div_tmp.dat, 0.0, u_out);
  op2_gemv(mesh, true, 1.0, DGConstants::DS, div_tmp.dat, 0.0, v_out);
  timer->endTimer("pr_proj - mult - drst");

  dg_dat_pool->releaseTempDatCells(div_tmp);

  timer->startTimer("pr_proj - mult - dir");
  op_par_loop(ins_2d_proj_5, "ins_2d_proj_5", mesh->cells,
              op_arg_dat(pen, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u_out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(v_out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("pr_proj - mult - dir");

  DGTempDat u_jump = dg_dat_pool->requestTempDatCells(DG_NPF * DG_NUM_FACES);
  DGTempDat v_jump = dg_dat_pool->requestTempDatCells(DG_NPF * DG_NUM_FACES);
  zero_dat(u_jump.dat);
  zero_dat(v_jump.dat);

  timer->startTimer("pr_proj - mult - indir");
  op_par_loop(ins_2d_proj_6, "ins_2d_proj_6", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(pen_f, -2, mesh->face2cells, 1, DG_FP_STR, OP_READ),
              op_arg_dat(u, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u_jump.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(v_jump.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  timer->endTimer("pr_proj - mult - indir");

  if(mesh->bface2cells) {
    op_par_loop(ins_2d_proj_7, "ins_2d_proj_7", mesh->bfaces,
                op_arg_gbl(&g0,   1, DG_FP_STR, OP_READ),
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(pen_f,   0, mesh->bface2cells, 1, DG_FP_STR, OP_READ),
                op_arg_dat(u, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u_jump.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(v_jump.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

  timer->startTimer("pr_proj - mult - emat");
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, u_jump.dat, 1.0, u_out);
  op2_gemv(mesh, false, 1.0, DGConstants::EMAT, v_jump.dat, 1.0, v_out);
  timer->endTimer("pr_proj - mult - emat");

  dg_dat_pool->releaseTempDatCells(u_jump);
  dg_dat_pool->releaseTempDatCells(v_jump);
  timer->endTimer("pr_proj - mult");
}

void INSSolverBase2D::project_velocity(op_dat dpdx, op_dat dpdy) {
  if(pr_projection_method == 1) {
    DGTempDat projRHS[2];
    projRHS[0] = dg_dat_pool->requestTempDatCells(DG_NP);
    projRHS[1] = dg_dat_pool->requestTempDatCells(DG_NP);
    // Calculate new velocity intermediate values
    op_par_loop(ins_2d_proj_rhs, "ins_2d_proj_rhs", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(projRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(projRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    mesh->mass(projRHS[0].dat);
    mesh->mass(projRHS[1].dat);

    DGTempDat proj_pen = dg_dat_pool->requestTempDatCells(1);
    DGTempDat proj_pen_f = dg_dat_pool->requestTempDatCells(1);
    DG_FP factor = dt * 1.0 * 1e3;
    // DG_FP factor = dt / Cr;
    // op_printf("Cr: %g\n", Cr);
    op_par_loop(ins_2d_proj_pen, "ins_2d_proj_pen", mesh->cells,
                op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen.dat, -1, OP_ID, 1, DG_FP_STR, OP_WRITE),
                op_arg_dat(proj_pen_f.dat, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));

    // Do the 2 vector linear solve with project_mat as the matrix
    int num_cells = 0;
    int num_converge = 0;
    DG_FP num_iter = 0.0;
    op_par_loop(ins_proj_cg_2d, "ins_proj_cg_2d", mesh->cells,
                op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen.dat, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(projRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(projRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_gbl(&num_cells, 1, "int", OP_INC),
                op_arg_gbl(&num_converge, 1, "int", OP_INC),
                op_arg_gbl(&num_iter, 1, DG_FP_STR, OP_INC));
    // op_printf("%d out of %d cells converged on projection step\n", num_converge, num_cells);
    // op_printf("Average iterations to converge on projection step %g\n", num_iter / (DG_FP)num_cells);
    if(num_cells != num_converge) {
      op_printf("%d out of %d cells converged on projection step\n", num_converge, num_cells);
      exit(-1);
    }

    dg_dat_pool->releaseTempDatCells(proj_pen);
    dg_dat_pool->releaseTempDatCells(proj_pen_f);
    dg_dat_pool->releaseTempDatCells(projRHS[0]);
    dg_dat_pool->releaseTempDatCells(projRHS[1]);
  } else if(pr_projection_method == 2) {
    DGTempDat projRHS[2];
    projRHS[0] = dg_dat_pool->requestTempDatCells(DG_NP);
    projRHS[1] = dg_dat_pool->requestTempDatCells(DG_NP);

    op_par_loop(ins_2d_proj_rhs, "ins_2d_proj_rhs", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(projRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(projRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    mesh->mass(projRHS[0].dat);
    mesh->mass(projRHS[1].dat);

    DGTempDat proj_pen = dg_dat_pool->requestTempDatCells(1);
    DGTempDat proj_pen_f = dg_dat_pool->requestTempDatCells(1);
    DG_FP factor = dt * 1.0 * 1e2;
    // DG_FP factor = dt / Cr;
    // op_printf("Cr: %g\n", Cr);
    op_par_loop(ins_2d_proj_pen, "ins_2d_proj_pen", mesh->cells,
                op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen.dat, -1, OP_ID, 1, DG_FP_STR, OP_WRITE),
                op_arg_dat(proj_pen_f.dat, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));

    DGTempDat cg_tmp[2];
    cg_tmp[0] = dg_dat_pool->requestTempDatCells(DG_NP);
    cg_tmp[1] = dg_dat_pool->requestTempDatCells(DG_NP);
    DGTempDat cg_r[2];
    cg_r[0] = dg_dat_pool->requestTempDatCells(DG_NP);
    cg_r[1] = dg_dat_pool->requestTempDatCells(DG_NP);
    DGTempDat cg_z[2];
    cg_z[0] = dg_dat_pool->requestTempDatCells(DG_NP);
    cg_z[1] = dg_dat_pool->requestTempDatCells(DG_NP);
    DGTempDat cg_p[2];
    cg_p[0] = dg_dat_pool->requestTempDatCells(DG_NP);
    cg_p[1] = dg_dat_pool->requestTempDatCells(DG_NP);

    int iter = 0;
    const int max_iter = 500;
    DG_FP residual = 0.0;
    const DG_FP tol = 1e-12;

    // Calculate first residual
    project_velocity_mat_mult(velTT[0], velTT[1], cg_tmp[0].dat,
            cg_tmp[1].dat, proj_pen.dat, proj_pen_f.dat);

    op_par_loop(ins_2d_proj_0, "ins_2d_proj_0", mesh->cells,
                op_arg_gbl(&residual, 1, DG_FP_STR, OP_INC),
                op_arg_dat(projRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(projRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    residual = sqrt(residual);
    // op_printf("residual: %g\n", residual);
    mesh->inv_mass(cg_z[0].dat);
    mesh->inv_mass(cg_z[1].dat);
    op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
                op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_p[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
                op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(cg_p[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    while(residual > tol && iter < max_iter) {
      project_velocity_mat_mult(cg_p[0].dat, cg_p[1].dat,
        cg_tmp[0].dat, cg_tmp[1].dat, proj_pen.dat, proj_pen_f.dat);

      DG_FP tmp_alpha_0 = 0.0;
      DG_FP tmp_alpha_1 = 0.0;
      op_par_loop(ins_2d_proj_1, "ins_2d_proj_1", mesh->cells,
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
      residual = 0.0;
      DG_FP tmp_beta_1 = 0.0;
      op_par_loop(ins_2d_proj_2, "ins_2d_proj_2", mesh->cells,
                  op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
                  op_arg_gbl(&residual, 1, DG_FP_STR, OP_INC),
                  op_arg_gbl(&tmp_beta_1, 1, DG_FP_STR, OP_INC),
                  op_arg_dat(cg_p[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_p[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

      residual = sqrt(residual);
      // op_printf("residual: %g\n", residual);

      if(residual < tol) break;

      mesh->inv_mass(cg_z[0].dat);
      mesh->inv_mass(cg_z[1].dat);

      DG_FP tmp_beta_0 = 0.0;
      op_par_loop(ins_2d_proj_3, "ins_2d_proj_3", mesh->cells,
                  op_arg_gbl(&tmp_beta_0, 1, DG_FP_STR, OP_INC),
                  op_arg_dat(cg_r[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_r[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

      DG_FP beta = tmp_beta_0 / tmp_beta_1;

      op_par_loop(ins_2d_proj_4, "ins_2d_proj_4", mesh->cells,
                  op_arg_gbl(&beta, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_z[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_z[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(cg_p[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(cg_p[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
      iter++;
    }

    // if(iter == max_iter)
    op_printf("residual: %g\n", residual);
    op_printf("iter: %d\n", iter);

    dg_dat_pool->releaseTempDatCells(cg_p[0]);
    dg_dat_pool->releaseTempDatCells(cg_p[1]);
    dg_dat_pool->releaseTempDatCells(cg_z[0]);
    dg_dat_pool->releaseTempDatCells(cg_z[1]);
    dg_dat_pool->releaseTempDatCells(cg_r[0]);
    dg_dat_pool->releaseTempDatCells(cg_r[1]);
    dg_dat_pool->releaseTempDatCells(cg_tmp[0]);
    dg_dat_pool->releaseTempDatCells(cg_tmp[1]);

    dg_dat_pool->releaseTempDatCells(proj_pen_f);
    dg_dat_pool->releaseTempDatCells(proj_pen);
    dg_dat_pool->releaseTempDatCells(projRHS[0]);
    dg_dat_pool->releaseTempDatCells(projRHS[1]);
  } else {
    // Calculate new velocity intermediate values
    op_par_loop(ins_pressure_update_2d, "ins_pressure_update_2d", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }
}

DG_FP INSSolverBase2D::max_vel() {
  DG_FP max_vel_tmp = 0.0;

  op_par_loop(ins_max_vel_2d, "ins_max_vel_2d", mesh->cells,
              op_arg_gbl(&max_vel_tmp, 1, DG_FP_STR, OP_MAX)
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  return std::max(1.0, sqrt(max_vel_tmp));
}

void INSSolverBase2D::add_to_pr_history() {
  const DG_FP t_n1 = time + dt;
  DGTempDat pr_copy = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
              op_arg_dat(pr, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pr_copy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  pr_history.push_back({t_n1, pr_copy});

  while(pr_history.size() > EXTRAPOLATE_HISTORY_SIZE) {
    dg_dat_pool->releaseTempDatCells(pr_history[0].second);
    pr_history.erase(pr_history.begin());
  }
}

DG_FP INSSolverBase2D::get_time() {
  return time;
}

DG_FP INSSolverBase2D::get_dt() {
  return dt;
}

void INSSolverBase2D::filter(op_dat in) {
  DGTempDat u_modal = dg_dat_pool->requestTempDatCells(DG_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, in, 0.0, u_modal.dat);
  op_par_loop(filter_2d, "filter_3d", mesh->cells,
              op_arg_gbl(&filter_alpha, 1, "double", OP_READ),
              op_arg_gbl(&filter_Nc, 1, "int", OP_READ),
              op_arg_gbl(&filter_sp,  1, "int", OP_READ),
              op_arg_dat(u_modal.dat, -1, OP_ID, DG_NP, "double", OP_RW));
  op2_gemv(mesh, false, 1.0, DGConstants::V, u_modal.dat, 0.0, in);
  dg_dat_pool->releaseTempDatCells(u_modal);
}

void INSSolverBase2D::calc_art_vis(op_dat in, op_dat out) {
  DGTempDat u_modal = dg_dat_pool->requestTempDatCells(DG_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, in, 0.0, u_modal.dat);

  op_par_loop(reset_tmp_node_dats, "reset_tmp_node_dats", mesh->nodes,
              op_arg_dat(nodes_data, -1, OP_ID, 1, DG_FP_STR, OP_WRITE),
              op_arg_dat(nodes_count, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));

  op_par_loop(modal_shock_detector_2d_0, "modal_shock_detector_2d_0", mesh->cells,
              op_arg_gbl(&shock_cap_max_diff, 1, DG_FP_STR, OP_READ),
              op_arg_dat(u_modal.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(nodes_data, -3, mesh->cell2nodes, 1, DG_FP_STR, OP_INC),
              op_arg_dat(nodes_count, -3, mesh->cell2nodes, 1, DG_FP_STR, OP_INC));

  dg_dat_pool->releaseTempDatCells(u_modal);

  op_par_loop(modal_shock_detector_2d_1, "modal_shock_detector_2d_1", mesh->cells,
              op_arg_dat(nodes_data, -3, mesh->cell2nodes, 1, DG_FP_STR, OP_READ),
              op_arg_dat(nodes_count, -3, mesh->cell2nodes, 1, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}

void INSSolverBase2D::dump_checkpoint_data(const std::string &filename) {
  op_fetch_data_hdf5_file(vel[0][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[0][1], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][1], filename.c_str());
  op_fetch_data_hdf5_file(n[0][0], filename.c_str());
  op_fetch_data_hdf5_file(n[0][1], filename.c_str());
  op_fetch_data_hdf5_file(n[1][0], filename.c_str());
  op_fetch_data_hdf5_file(n[1][1], filename.c_str());
  op_fetch_data_hdf5_file(dPdN[0], filename.c_str());
  op_fetch_data_hdf5_file(dPdN[1], filename.c_str());
  op_fetch_data_hdf5_file(pr, filename.c_str());

  op_write_const_hdf5("time", 1, "double", (char *)&time, filename.c_str());
  op_write_const_hdf5("prev_time", 1, "double", (char *)&prev_time, filename.c_str());
  op_write_const_hdf5("currentInd", 1, "int", (char *)&currentInd, filename.c_str());
  op_write_const_hdf5("it_pre_sub_cycle", 1, "int", (char *)&it_pre_sub_cycle, filename.c_str());
}

void INSSolverBase2D::dump_visualisation_data(const std::string &filename) {
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());

  if(values_to_save.count("velocity") != 0) {
    op_fetch_data_hdf5_file(vel[currentInd][0], filename.c_str());
    op_fetch_data_hdf5_file(vel[currentInd][1], filename.c_str());
  }

  if(values_to_save.count("non_linear") != 0) {
    op_fetch_data_hdf5_file(n[(currentInd + 1) % 2][0], filename.c_str());
    op_fetch_data_hdf5_file(n[(currentInd + 1) % 2][1], filename.c_str());
  }

  if(values_to_save.count("intermediate_velocities") != 0) {
    op_fetch_data_hdf5_file(velT[0], filename.c_str());
    op_fetch_data_hdf5_file(velT[1], filename.c_str());
    op_fetch_data_hdf5_file(velTT[0], filename.c_str());
    op_fetch_data_hdf5_file(velTT[1], filename.c_str());
  }

  if(values_to_save.count("pressure") != 0) {
    op_fetch_data_hdf5_file(pr, filename.c_str());
  }

  if(values_to_save.count("art_vis") != 0 && shock_capturing) {
    op_dat art_vis_x = op_decl_dat_temp(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "art_vis_x");
    calc_art_vis(velTT[0], art_vis_x);
    op_fetch_data_hdf5_file(art_vis_x, filename.c_str());
    op_free_dat_temp(art_vis_x);
    op_dat art_vis_y = op_decl_dat_temp(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "art_vis_y");
    calc_art_vis(velTT[1], art_vis_y);
    op_fetch_data_hdf5_file(art_vis_y, filename.c_str());
    op_free_dat_temp(art_vis_y);
  }
}

void INSSolverBase2D::zero_dat(op_dat dat) {
  if(dat->dim == DG_NP) {
    op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                op_arg_dat(dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  } else if(dat->dim == DG_NUM_FACES * DG_NPF) {
    op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
                op_arg_dat(dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  } else if(dat->dim == DG_NUM_FACES * DG_CUB_SURF_2D_NP) {
    op_par_loop(zero_cub_surf_2d, "zero_cub_surf_2d", mesh->cells,
                op_arg_dat(dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_WRITE));
  } else {
    dg_abort("Trying to zero dat with incompatible dimension");
  }
}

void INSSolverBase2D::update_time() {
  prev_time = time;
  time += dt;
}

// Getters
DGMesh2D* INSSolverBase2D::get_mesh() {
  return mesh;
}

op_dat INSSolverBase2D::get_vel_x() {
  return vel[currentInd][0];
}

op_dat INSSolverBase2D::get_vel_y() {
  return vel[currentInd][1];
}

op_dat INSSolverBase2D::get_pr() {
  return pr;
}

LinearSolver::Solvers INSSolverBase2D::set_solver_type(const std::string &str) {
  if(str == "petsc--amg") {
    return LinearSolver::PETSC_AMG;
  } else if(str == "jacobi") {
    return LinearSolver::PETSC_JACOBI;
  } else if(str == "block-jacobi") {
    return LinearSolver::PETSC_BLOCK_JACOBI;
  } else if(str == "inv-mass") {
    return LinearSolver::PETSC_INV_MASS;
  } else if(str == "p-multigrid") {
    return LinearSolver::PETSC_PMULTIGRID;
  } else {
    dg_abort("Unknown solver type: " + str);
    return LinearSolver::PETSC_AMG;
  }
}
