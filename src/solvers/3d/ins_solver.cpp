#include "solvers/3d/ins_solver.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"

extern DGConstants *constants;

#include "timing.h"
#include "config.h"
#include "linear_solvers/petsc_amg.h"
#include "linear_solvers/petsc_block_jacobi.h"
#include "linear_solvers/petsc_pmultigrid.h"
#include "linear_solvers/initial_guess_extrapolation.h"
#include "dg_dat_pool.h"
#include "dg_utils.h"

#include <string>
#include <iostream>
#include <stdexcept>

extern Timing *timer;
extern Config *config;
extern DGDatPool3D *dg_dat_pool;

#define ENSTROPY_FREQUENCY 10

INSSolver3D::INSSolver3D(DGMesh3D *m) {
  mesh = m;
  resuming = false;

  setup_common();

  std::string name;
  for(int i = 0; i < 3; i++) {
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

  dPdN[0] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_dPdN0");
  dPdN[1] = op_decl_dat(mesh->cells, 4 * DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_dPdN1");

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;

  dt = 0.0;
  time = 0.0;

  currentInd = 0;
}

INSSolver3D::INSSolver3D(DGMesh3D *m, const std::string &filename, const int iter) {
  mesh = m;
  resuming = true;

  setup_common();

  vel[0][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel00");
  vel[1][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel10");
  vel[0][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel01");
  vel[1][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel11");
  vel[0][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel02");
  vel[1][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_vel12");
  n[0][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n00");
  n[1][0] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n10");
  n[0][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n01");
  n[1][1] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n11");
  n[0][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n02");
  n[1][2] = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_n12");
  pr = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ins_solver_pr");
  dPdN[0] = op_decl_dat_hdf5(mesh->cells, 4 * DG_NPF, DG_FP_STR, filename.c_str(), "ins_solver_dPdN0");
  dPdN[1] = op_decl_dat_hdf5(mesh->cells, 4 * DG_NPF, DG_FP_STR, filename.c_str(), "ins_solver_dPdN1");

  currentInd = iter;

  if(iter > 0) {
    g0 = 1.5;
    a0 = 2.0;
    a1 = -0.5;
    b0 = 2.0;
    b1 = -1.0;
  } else {
    a0 = 1.0;
    a1 = 0.0;
    b0 = 1.0;
    b1 = 0.0;
    g0 = 1.0;
  }
}

void INSSolver3D::setup_common() {
  int tmp_div = 1;
  config->getInt("solver-options", "div_div", tmp_div);
  div_div_proj = tmp_div != 0;
  config->getInt("solver-options", "sub_cycle", sub_cycles);
  config->getInt("solver-options", "num_iter_before_sub_cycle", it_pre_sub_cycle);
  it_pre_sub_cycle = it_pre_sub_cycle > 1 ? it_pre_sub_cycle : 1;
  int tmp_eig = 1;
  config->getInt("solver-options", "extrapolate_initial_guess", tmp_eig);
  extrapolate_initial_guess = tmp_eig == 1;
  int tmp_shock = 1;
  config->getInt("solver-options", "shock_capturing", tmp_shock);
  shock_cap = tmp_shock == 1;

  std::string name;
  for(int i = 0; i < 3; i++) {
    name = "ins_solver_velT" + std::to_string(i);
    velT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "ins_solver_velTT" + std::to_string(i);
    velTT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
  }

  tmp_npf_bc = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_tmp_npf_bc");

  tmp_bc_1 = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_tmp_bc_1");
  bc_types = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_bc_types");

  if(div_div_proj) {
    proj_h   = op_decl_dat(mesh->cells, 1, DG_FP_STR, (DG_FP *)NULL, "ins_solver_proj_h");
  }

  pr_bc_types  = tmp_bc_1;
  vis_bc_types = tmp_bc_1;
  pr_bc  = tmp_npf_bc;
  vis_bc = tmp_npf_bc;

  pressureCoarseMatrix = new PoissonCoarseMatrix3D(mesh);
  pressureMatrix = new PoissonMatrixFreeDiag3D(mesh);
  viscosityMatrix = new MMPoissonMatrixFree3D(mesh);

  PETScPMultigrid *tmp_pressureSolver = new PETScPMultigrid(mesh);
  viscositySolver = new PETScInvMassSolver(mesh);

  int pr_tmp = 0;
  int vis_tmp = 0;
  config->getInt("solver-options", "pr_nullspace", pr_tmp);
  config->getInt("solver-options", "vis_nullspace", vis_tmp);
  tmp_pressureSolver->set_coarse_matrix(pressureCoarseMatrix);
  pressureSolver = tmp_pressureSolver;
  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_nullspace(pr_tmp == 1);
  viscositySolver->set_matrix(viscosityMatrix);
  viscositySolver->set_nullspace(vis_tmp == 1);
}

INSSolver3D::~INSSolver3D() {
  delete pressureCoarseMatrix;
  delete pressureMatrix;
  delete viscosityMatrix;
  delete (PETScPMultigrid *)pressureSolver;
  delete viscositySolver;
}

void INSSolver3D::init(const DG_FP re, const DG_FP refVel) {
  timer->startTimer("INSSolver3D - Init");
  reynolds = re;

  // Set initial conditions
  if(!resuming) {
    op_par_loop(ins_3d_set_ic, "ins_3d_set_ic", mesh->cells,
                op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[0][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }

  h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  op_printf("h: %g\n", h);

  double tmp_dt = -1.0;
  config->getDouble("solver-options", "dt", tmp_dt);
  forced_dt = tmp_dt > 0.0;
  if(forced_dt) {
    dt = tmp_dt;
  } else {
    sub_cycle_dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1) * max_vel());
    dt = sub_cycle_dt;
    if(resuming)
      dt = sub_cycles > 1 ? sub_cycle_dt * sub_cycles : sub_cycle_dt;
  }
  // dt *= 1e-2;
  op_printf("INS dt is %g\n", dt);
  time = dt * currentInd;
  currentInd = currentInd % 2;

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_bc_types, "ins_3d_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -3, mesh->bface2nodes, 3, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }

  // Set up pressure projection
  if(div_div_proj) {
    DGTempDat tmp_npf = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
    op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
                op_arg_dat(tmp_npf.dat, -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));

    op_par_loop(ins_3d_proj_setup_0, "ins_3d_proj_setup_0", mesh->faces,
                op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
                op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_npf.dat, -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));

    op_par_loop(ins_3d_proj_setup_1, "ins_3d_proj_setup_1", mesh->cells,
                op_arg_dat(tmp_npf.dat, -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_READ),
                op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));
    dg_dat_pool->releaseTempDatCells(tmp_npf);
  }

  pressureSolver->init();
  viscositySolver->init();

  enstropy_counter = ENSTROPY_FREQUENCY;
  record_enstrophy();

  timer->endTimer("INSSolver3D - Init");
}

void INSSolver3D::step() {
  timer->startTimer("INSSolver3D - Advection");
  advection();
  timer->endTimer("INSSolver3D - Advection");

  timer->startTimer("INSSolver3D - Shock Capturing");
  if(shock_cap) {
    shock_capture_filter_dat(velT[0]);
    shock_capture_filter_dat(velT[1]);
    shock_capture_filter_dat(velT[2]);
  }
  timer->endTimer("INSSolver3D - Shock Capturing");

  timer->startTimer("INSSolver3D - Pressure");
  pressure();
  timer->endTimer("INSSolver3D - Pressure");

  timer->startTimer("INSSolver3D - Viscosity");
  viscosity();
  timer->endTimer("INSSolver3D - Viscosity");

  currentInd = (currentInd + 1) % 2;
  time += dt;
  g0 = 1.5;
  a0 = 2.0;
  a1 = -0.5;
  b0 = 2.0;
  b1 = -1.0;

  if(!forced_dt) {
    if(it_pre_sub_cycle > 1) {
      it_pre_sub_cycle--;
    } else {
      sub_cycle_dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1) * max_vel());
      dt = sub_cycles > 1 ? sub_cycle_dt * sub_cycles : sub_cycle_dt;
      it_pre_sub_cycle = 0;
    }
  }

  record_enstrophy();
}

DG_FP INSSolver3D::max_vel() {
  DG_FP max_vel_tmp = 0.0;

  op_par_loop(ins_3d_max_vel, "ins_3d_max_vel", mesh->cells,
              op_arg_gbl(&max_vel_tmp, 1, DG_FP_STR, OP_MAX)
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  return std::max(1.0, sqrt(max_vel_tmp));
}

void INSSolver3D::advection() {
  if(time == 0.0 || sub_cycles < 1 || it_pre_sub_cycle != 0) {
    advec_standard();
  } else {
    advec_sub_cycle();
  }
}

void INSSolver3D::advec_current_non_linear() {
  DGTempDat f[3][3];
  f[0][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[0][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[0][2] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][2] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[2][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[2][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[2][2] = dg_dat_pool->requestTempDatCells(DG_NP);
  op_par_loop(ins_3d_advec_0, "ins_3d_advec_0", mesh->cells,
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div(f[0][0].dat, f[0][1].dat, f[0][2].dat, n[currentInd][0]);
  mesh->div(f[1][0].dat, f[1][1].dat, f[1][2].dat, n[currentInd][1]);
  mesh->div(f[2][0].dat, f[2][1].dat, f[2][2].dat, n[currentInd][2]);

  dg_dat_pool->releaseTempDatCells(f[0][0]);
  dg_dat_pool->releaseTempDatCells(f[0][1]);
  dg_dat_pool->releaseTempDatCells(f[0][2]);
  dg_dat_pool->releaseTempDatCells(f[1][0]);
  dg_dat_pool->releaseTempDatCells(f[1][1]);
  dg_dat_pool->releaseTempDatCells(f[1][2]);
  dg_dat_pool->releaseTempDatCells(f[2][0]);
  dg_dat_pool->releaseTempDatCells(f[2][1]);
  dg_dat_pool->releaseTempDatCells(f[2][2]);

  DGTempDat tmp_advec_flux0 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_advec_flux1 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_advec_flux2 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(tmp_advec_flux0.dat, -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux1.dat, -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux2.dat, -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));

  // Flux across faces
  op_par_loop(ins_3d_advec_1, "ins_3d_advec_1", mesh->faces,
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_advec_flux0.dat, -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(tmp_advec_flux1.dat, -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(tmp_advec_flux2.dat, -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));

  // Boundary flux
  if(mesh->bface2cells) {
    op_par_loop(ins_3d_advec_2, "ins_3d_advec_2", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_advec_flux0.dat, 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_advec_flux1.dat, 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_advec_flux2.dat, 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux0.dat, 1.0, n[currentInd][0]);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux1.dat, 1.0, n[currentInd][1]);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux2.dat, 1.0, n[currentInd][2]);

  dg_dat_pool->releaseTempDatCells(tmp_advec_flux0);
  dg_dat_pool->releaseTempDatCells(tmp_advec_flux1);
  dg_dat_pool->releaseTempDatCells(tmp_advec_flux2);
}

void INSSolver3D::advec_standard() {
  advec_current_non_linear();

  // Calculate the intermediate velocity values
  op_par_loop(ins_3d_advec_3, "ins_3d_advec_3", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(n[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}

void INSSolver3D::advec_sub_cycle_rk_step(const DG_FP time_sc, op_dat u, op_dat v, op_dat w) {
  // Request temporary dats
  DGTempDat advec_sc_rk[3][3];
  advec_sc_rk[0][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[0][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[0][2] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[1][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[1][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[1][2] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[2][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[2][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_rk[2][2] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_3d_advec_sc_copy, "ins_3d_advec_sc_copy", mesh->cells,
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(w, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc_rk[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc_rk[0][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  for(int rk_step = 0; rk_step < 3; rk_step++) {
    double rk_time = time_sc;
    if(rk_step == 1) rk_time += sub_cycle_dt;
    if(rk_step == 2) rk_time += 0.5 * sub_cycle_dt;
    const int rk_ind = rk_step == 2 ? 2 : rk_step + 1;
    advec_sub_cycle_rhs(advec_sc_rk[0][0].dat, advec_sc_rk[0][1].dat, advec_sc_rk[0][2].dat,
                        advec_sc_rk[rk_ind][0].dat, advec_sc_rk[rk_ind][1].dat,
                        advec_sc_rk[rk_ind][2].dat, rk_time);
    // Set up next step
    if(rk_step == 0) {
      op_par_loop(ins_3d_advec_sc_rk_0, "ins_3d_advec_sc_rk_0", mesh->cells,
                  op_arg_gbl(&sub_cycle_dt, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(w, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                  op_arg_dat(advec_sc_rk[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                  op_arg_dat(advec_sc_rk[0][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    } else if(rk_step == 1) {
      op_par_loop(ins_3d_advec_sc_rk_1, "ins_3d_advec_sc_rk_1", mesh->cells,
                  op_arg_gbl(&sub_cycle_dt, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(w, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(advec_sc_rk[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(advec_sc_rk[1][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                  op_arg_dat(advec_sc_rk[2][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[2][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[2][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(advec_sc_rk[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                  op_arg_dat(advec_sc_rk[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                  op_arg_dat(advec_sc_rk[0][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    }
  }
  // Update velT
  op_par_loop(ins_3d_advec_sc_rk_2, "ins_3d_advec_sc_rk_2", mesh->cells,
              op_arg_gbl(&sub_cycle_dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[1][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[2][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[2][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_rk[2][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(v, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(w, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  dg_dat_pool->releaseTempDatCells(advec_sc_rk[0][0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[0][1]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[0][2]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[1][0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[1][1]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[1][2]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[2][0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[2][1]);
  dg_dat_pool->releaseTempDatCells(advec_sc_rk[2][2]);
}

void INSSolver3D::advec_sub_cycle() {
  op_par_loop(ins_3d_advec_sc_copy, "ins_3d_advec_sc_copy", mesh->cells,
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(velT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  // Advance 2 * number of subcycles
  for(int i = 0; i < 2 * sub_cycles; i++) {
    advec_sub_cycle_rk_step(time - dt + i * sub_cycle_dt, velT[0], velT[1], velT[2]);
  }

  DGTempDat advec_sc_tmp[3];
  advec_sc_tmp[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_tmp[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc_tmp[2] = dg_dat_pool->requestTempDatCells(DG_NP);

  // Other velocity
  op_par_loop(ins_3d_advec_sc_copy, "ins_3d_advec_sc_copy", mesh->cells,
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc_tmp[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  // Advance number of subcycles
  for(int i = 0; i < sub_cycles; i++) {
    advec_sub_cycle_rk_step(time + i * sub_cycle_dt, advec_sc_tmp[0].dat, advec_sc_tmp[1].dat, advec_sc_tmp[2].dat);
  }

  // Get final velT
  op_par_loop(ins_3d_advec_sc_update, "ins_3d_advec_sc_update", mesh->cells,
              op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_tmp[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_tmp[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc_tmp[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(velT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  dg_dat_pool->releaseTempDatCells(advec_sc_tmp[0]);
  dg_dat_pool->releaseTempDatCells(advec_sc_tmp[1]);
  dg_dat_pool->releaseTempDatCells(advec_sc_tmp[2]);

  // Needed for Pressure boundary conditions
  advec_current_non_linear();
}

void INSSolver3D::advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat w_in,
                                      op_dat u_out, op_dat v_out, op_dat w_out,
                                      const double t) {
  double t0 = time - dt;
  double t1 = time;
  double tI = t;
  DGTempDat advec_sc[3];
  advec_sc[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  advec_sc[2] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_3d_advec_sc_interp, "ins_3d_advec_sc_interp", mesh->cells,
              op_arg_gbl(&t0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&t1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&tI, 1, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(advec_sc[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  DGTempDat f[3][3];
  f[0][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[0][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[0][2] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[1][2] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[2][0] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[2][1] = dg_dat_pool->requestTempDatCells(DG_NP);
  f[2][2] = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(ins_3d_advec_sc_rhs_0, "ins_3d_advec_sc_rhs_0", mesh->cells,
              op_arg_dat(u_in, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v_in, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(w_in, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f[0][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[0][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[1][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(f[2][2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div(f[0][0].dat, f[0][1].dat, f[0][2].dat, u_out);
  mesh->div(f[1][0].dat, f[1][1].dat, f[1][2].dat, v_out);
  mesh->div(f[2][0].dat, f[2][1].dat, f[2][2].dat, w_out);

  dg_dat_pool->releaseTempDatCells(f[0][0]);
  dg_dat_pool->releaseTempDatCells(f[0][1]);
  dg_dat_pool->releaseTempDatCells(f[0][2]);
  dg_dat_pool->releaseTempDatCells(f[1][0]);
  dg_dat_pool->releaseTempDatCells(f[1][1]);
  dg_dat_pool->releaseTempDatCells(f[1][2]);
  dg_dat_pool->releaseTempDatCells(f[2][0]);
  dg_dat_pool->releaseTempDatCells(f[2][1]);
  dg_dat_pool->releaseTempDatCells(f[2][2]);

  DGTempDat tmp_advec_flux0 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_advec_flux1 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat tmp_advec_flux2 = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_3, "zero_npf_3", mesh->cells,
              op_arg_dat(tmp_advec_flux0.dat, -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux1.dat, -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux2.dat, -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(ins_3d_advec_sc_rhs_1, "ins_3d_advec_sc_rhs_1", mesh->faces,
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(u_in, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v_in, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(w_in, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[0].dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[1].dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(advec_sc[2].dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_advec_flux0.dat, -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(tmp_advec_flux1.dat, -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
              op_arg_dat(tmp_advec_flux2.dat, -2, mesh->face2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_advec_sc_rhs_2, "ins_3d_advec_sc_rhs_2", mesh->bfaces,
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(w_in, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc[0].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc[1].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(advec_sc[2].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_advec_flux0.dat, 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_advec_flux1.dat, 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_advec_flux2.dat, 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));
  }

  dg_dat_pool->releaseTempDatCells(advec_sc[0]);
  dg_dat_pool->releaseTempDatCells(advec_sc[1]);
  dg_dat_pool->releaseTempDatCells(advec_sc[2]);

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux0.dat, 1.0, u_out);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux1.dat, 1.0, v_out);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_advec_flux2.dat, 1.0, w_out);

  dg_dat_pool->releaseTempDatCells(tmp_advec_flux0);
  dg_dat_pool->releaseTempDatCells(tmp_advec_flux1);
  dg_dat_pool->releaseTempDatCells(tmp_advec_flux2);
}

void INSSolver3D::pressure() {
  timer->startTimer("INSSolver3D - Pressure RHS");
  DGTempDat divVelT = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat curlVel[3];
  curlVel[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  curlVel[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  curlVel[2] = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat curl2Vel[3];
  curl2Vel[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  curl2Vel[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  curl2Vel[2] = dg_dat_pool->requestTempDatCells(DG_NP);
  // mesh->div(velT[0], velT[1], velT[2], divVelT);
  mesh->div_with_central_flux(velT[0], velT[1], velT[2], divVelT.dat);
  mesh->curl(vel[currentInd][0], vel[currentInd][1], vel[currentInd][2],
             curlVel[0].dat, curlVel[1].dat, curlVel[2].dat);
  mesh->curl(curlVel[0].dat, curlVel[1].dat, curlVel[2].dat, curl2Vel[0].dat,
             curl2Vel[1].dat, curl2Vel[2].dat);

  dg_dat_pool->releaseTempDatCells(curlVel[0]);
  dg_dat_pool->releaseTempDatCells(curlVel[1]);
  dg_dat_pool->releaseTempDatCells(curlVel[2]);

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_pr_0, "ins_3d_pr_0", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(curl2Vel[0].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(curl2Vel[1].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(curl2Vel[2].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dPdN[currentInd], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));
  }

  dg_dat_pool->releaseTempDatCells(curl2Vel[0]);
  dg_dat_pool->releaseTempDatCells(curl2Vel[1]);
  dg_dat_pool->releaseTempDatCells(curl2Vel[2]);

  op_par_loop(ins_3d_pr_1, "ins_3d_pr_1", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(divVelT.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT.dat);
  mesh->mass(divVelT.dat);

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_pr_2, "ins_3d_pr_2", mesh->bfaces,
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(pr_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(pr_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }
  timer->endTimer("INSSolver3D - Pressure RHS");

  timer->startTimer("INSSolver3D - Pressure Linear Solve");
  pressureCoarseMatrix->set_bc_types(pr_bc_types);
  pressureMatrix->set_bc_types(pr_bc_types);
  pressureSolver->set_bcs(pr_bc);

  if(extrapolate_initial_guess)
    initial_guess_extrapolation(mesh, pr_history, pr, time + dt);

  bool converged = pressureSolver->solve(divVelT.dat, pr);
  if(!converged)
    throw std::runtime_error("\nPressure solve failed to converge\n");

  if(extrapolate_initial_guess)
    add_to_pr_history();
  timer->endTimer("INSSolver3D - Pressure Linear Solve");

  dg_dat_pool->releaseTempDatCells(divVelT);

  timer->startTimer("INSSolver3D - Pressure Projection");
  project_velocity();
  timer->endTimer("INSSolver3D - Pressure Projection");
}

void INSSolver3D::project_velocity() {
  DGTempDat dpdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdy = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdz = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->grad_with_central_flux(pr, dpdx.dat, dpdy.dat, dpdz.dat);

  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));

  if(div_div_proj) {
    DGTempDat projRHS[3];
    projRHS[0] = dg_dat_pool->requestTempDatCells(DG_NP);
    projRHS[1] = dg_dat_pool->requestTempDatCells(DG_NP);
    projRHS[2] = dg_dat_pool->requestTempDatCells(DG_NP);

    op_par_loop(ins_3d_proj_rhs, "ins_3d_proj_rhs", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdz.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(projRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(projRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(projRHS[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    mesh->mass(projRHS[0].dat);
    mesh->mass(projRHS[1].dat);
    mesh->mass(projRHS[2].dat);

    DGTempDat proj_pen = dg_dat_pool->requestTempDatCells(1);
    DG_FP factor = dt * 1.0;
    // DG_FP factor = dt / Cr;
    // op_printf("Cr: %g\n", Cr);
    op_par_loop(ins_3d_proj_pen, "ins_3d_proj_pen", mesh->cells,
                op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a0, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&a1, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, -1, OP_ID, 10, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[(currentInd + 1) % 2][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen.dat, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));

    int num_cells = 0;
    int num_converge = 0;
    DG_FP num_iter = 0.0;
    op_par_loop(ins_3d_proj_cg, "ins_3d_proj_cg", mesh->cells,
                op_arg_dat(mesh->rx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->tx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->ry, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sy, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->ty, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->rz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->sz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->tz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->J,  -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen.dat, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(projRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(projRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(projRHS[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(velTT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_gbl(&num_cells, 1, "int", OP_INC),
                op_arg_gbl(&num_converge, 1, "int", OP_INC),
                op_arg_gbl(&num_iter, 1, DG_FP_STR, OP_INC));

    dg_dat_pool->releaseTempDatCells(proj_pen);
    dg_dat_pool->releaseTempDatCells(projRHS[0]);
    dg_dat_pool->releaseTempDatCells(projRHS[1]);
    dg_dat_pool->releaseTempDatCells(projRHS[2]);
    if(num_cells != num_converge) {
      op_printf("%d out of %d cells converged on projection step\n", num_converge, num_cells);
      exit(-1);
    }
  } else {
    op_par_loop(ins_3d_pr_3, "ins_3d_pr_3", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdz.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }

  dg_dat_pool->releaseTempDatCells(dpdx);
  dg_dat_pool->releaseTempDatCells(dpdy);
  dg_dat_pool->releaseTempDatCells(dpdz);
}

void INSSolver3D::viscosity() {
  timer->startTimer("INSSolver3D - Viscosity RHS");
  DGTempDat visRHS[3];
  visRHS[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  visRHS[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  visRHS[2] = dg_dat_pool->requestTempDatCells(DG_NP);

  DG_FP factor = reynolds / dt;
  op_par_loop(ins_3d_vis_0, "ins_3d_vis_0", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visRHS[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->mass(visRHS[0].dat);
  mesh->mass(visRHS[1].dat);
  mesh->mass(visRHS[2].dat);
  timer->endTimer("INSSolver3D - Viscosity RHS");

  DG_FP vis_time = time + dt;
  if(mesh->bface2cells) {
    op_par_loop(ins_3d_vis_x, "ins_3d_vis_x", mesh->bfaces,
                op_arg_gbl(&vis_time, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&g0, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  timer->startTimer("INSSolver3D - Viscosity Linear Solve");
  factor = g0 * reynolds / dt;
  if(factor != viscosityMatrix->get_factor()) {
    viscosityMatrix->set_factor(factor);
    viscosityMatrix->set_bc_types(vis_bc_types);
    // viscosityMatrix->calc_mat();
    viscositySolver->setFactor(1.0 / factor);
  }
  viscositySolver->set_bcs(vis_bc);
  bool convergedX = viscositySolver->solve(visRHS[0].dat, vel[(currentInd + 1) % 2][0]);
  if(!convergedX)
    throw std::runtime_error("\nViscosity X solve failed to converge\n");

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_vis_y, "ins_3d_vis_y", mesh->bfaces,
                op_arg_gbl(&vis_time, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&g0, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  viscositySolver->set_bcs(vis_bc);
  bool convergedY = viscositySolver->solve(visRHS[1].dat, vel[(currentInd + 1) % 2][1]);
  if(!convergedY)
    throw std::runtime_error("\nViscosity Y solve failed to converge\n");

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_vis_z, "ins_3d_vis_z", mesh->bfaces,
                op_arg_gbl(&vis_time, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&g0, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[2], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  viscositySolver->set_bcs(vis_bc);
  bool convergedZ = viscositySolver->solve(visRHS[2].dat, vel[(currentInd + 1) % 2][2]);
  if(!convergedZ)
    throw std::runtime_error("\nViscosity Z solve failed to converge\n");

  dg_dat_pool->releaseTempDatCells(visRHS[0]);
  dg_dat_pool->releaseTempDatCells(visRHS[1]);
  dg_dat_pool->releaseTempDatCells(visRHS[2]);
  timer->endTimer("INSSolver3D - Viscosity Linear Solve");
}

DG_FP INSSolver3D::get_time() {
  return time;
}

DG_FP INSSolver3D::get_dt() {
  return dt;
}

void INSSolver3D::dump_data(const std::string &filename) {
  timer->startTimer("INSSolver3D - Dump Data");
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
  op_fetch_data_hdf5_file(mesh->z, filename.c_str());
  op_fetch_data_hdf5_file(vel[0][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[0][1], filename.c_str());
  op_fetch_data_hdf5_file(vel[0][2], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][0], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][1], filename.c_str());
  op_fetch_data_hdf5_file(vel[1][2], filename.c_str());
  op_fetch_data_hdf5_file(n[0][0], filename.c_str());
  op_fetch_data_hdf5_file(n[0][1], filename.c_str());
  op_fetch_data_hdf5_file(n[0][2], filename.c_str());
  op_fetch_data_hdf5_file(n[1][0], filename.c_str());
  op_fetch_data_hdf5_file(n[1][1], filename.c_str());
  op_fetch_data_hdf5_file(n[1][2], filename.c_str());
  op_fetch_data_hdf5_file(velT[0], filename.c_str());
  op_fetch_data_hdf5_file(velT[1], filename.c_str());
  op_fetch_data_hdf5_file(velT[2], filename.c_str());
  op_fetch_data_hdf5_file(velTT[0], filename.c_str());
  op_fetch_data_hdf5_file(velTT[1], filename.c_str());
  op_fetch_data_hdf5_file(velTT[2], filename.c_str());
  op_fetch_data_hdf5_file(pr, filename.c_str());
  op_fetch_data_hdf5_file(dPdN[0], filename.c_str());
  op_fetch_data_hdf5_file(dPdN[1], filename.c_str());
  timer->endTimer("INSSolver3D - Dump Data");
}

void INSSolver3D::add_to_pr_history() {
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

DG_FP INSSolver3D::calc_enstrophy() {
  DGTempDat curl[3];
  curl[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  curl[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  curl[2] = dg_dat_pool->requestTempDatCells(DG_NP);

  mesh->curl(vel[currentInd][0], vel[currentInd][1], vel[currentInd][2],
             curl[0].dat, curl[1].dat, curl[2].dat);

  op_par_loop(ins_3d_enstrophy_0, "ins_3d_enstrophy_0", mesh->cells,
              op_arg_dat(curl[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(curl[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(curl[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  mesh->mass(curl[2].dat);

  DG_FP enstropy = 0.0;

  op_par_loop(ins_3d_enstrophy_1, "ins_3d_enstrophy_1", mesh->cells,
              op_arg_gbl(&enstropy, 1, DG_FP_STR, OP_INC),
              op_arg_dat(curl[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(curl[0]);
  dg_dat_pool->releaseTempDatCells(curl[1]);
  dg_dat_pool->releaseTempDatCells(curl[2]);

  // Multiply by mu and divide by volume
  return 0.000625 * enstropy / 31.006276680299820175476315;
}

void INSSolver3D::record_enstrophy() {
  if(enstropy_counter % ENSTROPY_FREQUENCY == 0) {
    enstropy_history.push_back({time, calc_enstrophy()});
  }
  enstropy_counter++;
}

#include <fstream>
#include <iostream>

#include <iomanip>
#include <sstream>

std::string doubleToText(const double &d) {
    std::stringstream ss;
    ss << std::setprecision(15);
    ss << d;
    return ss.str();
}

#ifdef INS_MPI
#include "mpi.h"
#endif

void INSSolver3D::save_enstropy_history(const std::string &filename) {
  #ifdef INS_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
  #endif
  std::ofstream file(filename);

  file << "time,enstropy" << std::endl;

  for(int i = 0; i < enstropy_history.size(); i++) {
    file << doubleToText(enstropy_history[i].first) << ",";
    file << doubleToText(enstropy_history[i].second) << std::endl;
  }

  file.close();
  #ifdef INS_MPI
  }
  #endif
}

void INSSolver3D::shock_capture_filter_dat(op_dat in) {
  DGTempDat u_hat = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat u_modal = dg_dat_pool->requestTempDatCells(DG_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, in, 0.0, u_modal.dat);

  const double *r_ptr = constants->get_mat_ptr(DGConstants::R) + (DG_ORDER - 1) * DG_NP;
  const double *s_ptr = constants->get_mat_ptr(DGConstants::S) + (DG_ORDER - 1) * DG_NP;
  const double *t_ptr = constants->get_mat_ptr(DGConstants::T) + (DG_ORDER - 1) * DG_NP;

  std::vector<DG_FP> r_vec, s_vec, t_vec;
  for(int i = 0; i < DG_NP; i++) {
    r_vec.push_back(r_ptr[i]);
    s_vec.push_back(s_ptr[i]);
    t_vec.push_back(t_ptr[i]);
  }

  std::vector<DG_FP> simplex_vals = DGUtils::val_at_pt_N_1_3d_get_simplexes(r_vec, s_vec, t_vec, DG_ORDER);

  op_par_loop(discont_sensor_0, "discont_sensor_0", mesh->cells,
              op_arg_gbl(simplex_vals.data(), DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(in, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_modal.dat, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_hat.dat, -1, OP_ID, DG_NP, "double", OP_WRITE));

  double max_alpha = 18;
  // double e0 = h;
  // double s0 = 1.0 / (double)(DG_ORDER * DG_ORDER * DG_ORDER * DG_ORDER);
  double s0 = -2.0;
  // double k  = 5.0;
  double k = 1.0;
  op_par_loop(discont_sensor_filter, "discont_sensor_filter", mesh->cells,
              op_arg_gbl(&max_alpha, 1, "double", OP_READ),
              op_arg_gbl(&s0, 1, "double", OP_READ),
              op_arg_gbl(&k,  1, "double", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 10, "double", OP_READ),
              op_arg_dat(in, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_hat.dat, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u_modal.dat, -1, OP_ID, DG_NP, "double", OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::V, u_modal.dat, 0.0, in);

  dg_dat_pool->releaseTempDatCells(u_hat);
  dg_dat_pool->releaseTempDatCells(u_modal);
}
