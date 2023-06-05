#include "solvers/2d/ins_solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"

#include "timing.h"
#include "config.h"
#include "linear_solvers/petsc_amg.h"
#include "linear_solvers/petsc_block_jacobi.h"
#include "linear_solvers/petsc_pmultigrid.h"

extern Timing *timer;
extern Config *config;
extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

using namespace std;

INSSolver2D::INSSolver2D(DGMesh2D *m) {
  mesh = m;
  resuming = false;

  setup_common();

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

  currentInd = 0;

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;
}

INSSolver2D::INSSolver2D(DGMesh2D *m, const std::string &filename, const int iter) {
  mesh = m;
  resuming = true;

  setup_common();

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

void INSSolver2D::setup_common() {
  pressureMatrix = new PoissonMatrixFreeDiag2D(mesh);
  pressureCoarseMatrix = new PoissonCoarseMatrix2D(mesh);
  // viscosityMatrix = new MMPoissonMatrixOverInt2D(mesh);
  // viscosityMatrix = new MMPoissonMatrixFreeOverInt2D(mesh);
  PETScPMultigrid *tmp_pressureSolver = new PETScPMultigrid(mesh);
  tmp_pressureSolver->set_coarse_matrix(pressureCoarseMatrix);
  pressureSolver = tmp_pressureSolver;
  // viscositySolver = new PETScBlockJacobiSolver(mesh);
  // viscositySolver = new PETScInvMassSolver(mesh);

  int pr_tmp = 0;
  int vis_tmp = 0;
  config->getInt("solver-options", "pr_nullspace", pr_tmp);
  config->getInt("solver-options", "vis_nullspace", vis_tmp);

  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_nullspace(pr_tmp == 1);
  // viscositySolver->set_matrix(viscosityMatrix);
  // viscositySolver->set_nullspace(vis_tmp == 1);

  std::string name;
  for(int i = 0; i < 2; i++) {
    name = "ins_solver_velT" + std::to_string(i);
    velT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
    name = "ins_solver_velTT" + std::to_string(i);
    velTT[i] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, name.c_str());
  }

  bc_types     = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_bc_types");
  pr_bc_types  = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_pr_bc_types");
  vis_bc_types = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_vis_bc_types");

  // proj_pen = op_decl_dat(mesh->cells, 1, DG_FP_STR, (DG_FP *)NULL, "proj_pen");
  // proj_h   = op_decl_dat(mesh->cells, 1, DG_FP_STR, (DG_FP *)NULL, "proj_h");
}

INSSolver2D::~INSSolver2D() {
  delete pressureCoarseMatrix;
  delete pressureMatrix;
  // delete viscosityMatrix;
  delete pressureSolver;
  // delete viscositySolver;
}

void INSSolver2D::init(const DG_FP re, const DG_FP refVel) {
  timer->startTimer("INSSolver2D - Init");
  reynolds = re;

  // Set initial conditions
  if(!resuming) {
    op_par_loop(ins_2d_set_ic, "ins_2d_set_ic", mesh->cells,
                op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    op_par_loop(zero_npf_2, "zero_npf_2", mesh->cells,
                op_arg_dat(dPdN[0], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
                op_arg_dat(dPdN[1], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  }

  dt = numeric_limits<DG_FP>::max();
  op_par_loop(calc_dt, "calc_dt", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_MIN));
  dt = dt / (DG_ORDER * DG_ORDER * refVel);
  op_printf("dt: %g\n", dt);

  time = dt * currentInd;
  currentInd = currentInd % 2;

  if(mesh->bface2nodes) {
    op_par_loop(ins_bc_types, "ins_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -3, mesh->bface2nodes, 3, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,     -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(pr_bc_types,  -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }

  pressureCoarseMatrix->set_bc_types(pr_bc_types);
  pressureMatrix->set_bc_types(pr_bc_types);
  pressureSolver->init();
  // viscositySolver->init();

  timer->endTimer("INSSolver2D - Init");
}

void INSSolver2D::step() {
  timer->startTimer("INSSolver2D - Advection");
  advection();
  timer->endTimer("INSSolver2D - Advection");

  timer->startTimer("INSSolver2D - Pressure");
  pressure();
  timer->endTimer("INSSolver2D - Pressure");

  // timer->startTimer("Shock Capturing");
  // shock_capturing();
  // timer->endTimer("Shock Capturing");

  // timer->startTimer("INSSolver2D - Viscosity");
  // viscosity();
  // timer->endTimer("INSSolver2D - Viscosity");

  currentInd = (currentInd + 1) % 2;
  time += dt;
  g0 = 1.5;
  a0 = 2.0;
  a1 = -0.5;
  b0 = 2.0;
  b1 = -1.0;
}

// Calculate Nonlinear Terms
void INSSolver2D::advection() {
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

  op_par_loop(zero_npf_2, "zero_npf_2", mesh->cells,
              op_arg_dat(tmp_advec_flux0.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_advec_flux1.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

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

  // Calculate the intermediate velocity values
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

bool INSSolver2D::pressure() {
  timer->startTimer("INSSolver2D - Pressure RHS");
  DGTempDat divVelT = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat curlVel = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat gradCurlVel[2];
  gradCurlVel[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  gradCurlVel[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->div_with_central_flux(velT[0], velT[1], divVelT.dat);
  mesh->curl(vel[currentInd][0], vel[currentInd][1], curlVel.dat);
  mesh->grad(curlVel.dat, gradCurlVel[0].dat, gradCurlVel[1].dat);

  dg_dat_pool->releaseTempDatCells(curlVel);

  // Apply Neumann pressure boundary conditions
  if(mesh->bface2cells) {
    op_par_loop(ins_pressure_bc_2d, "ins_pressure_bc_2d", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(n[currentInd][1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gradCurlVel[0].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(gradCurlVel[1].dat, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dPdN[currentInd],   0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

  dg_dat_pool->releaseTempDatCells(gradCurlVel[0]);
  dg_dat_pool->releaseTempDatCells(gradCurlVel[1]);

  DGTempDat prBC = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  // Apply Dirichlet BCs
  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(prBC.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  // Calculate RHS of pressure solve
  // This assumes that the boundaries will always be order DG_ORDER
  op_par_loop(ins_pressure_rhs_2d, "ins_pressure_rhs_2d", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_RW),
              op_arg_dat(divVelT.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  mesh->mass(divVelT.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT.dat);
  timer->endTimer("INSSolver2D - Pressure RHS");

  // Call PETSc linear solver
  timer->startTimer("INSSolver2D - Pressure Linear Solve");
  pressureSolver->set_bcs(prBC.dat);
  bool converged = pressureSolver->solve(divVelT.dat, pr);
  dg_dat_pool->releaseTempDatCells(divVelT);
  dg_dat_pool->releaseTempDatCells(prBC);
  timer->endTimer("INSSolver2D - Pressure Linear Solve");

  timer->startTimer("INSSolver2D - Pressure Projection");
  project_velocity();
  timer->endTimer("INSSolver2D - Pressure Projection");

  return converged;
}

void INSSolver2D::project_velocity() {
  /*
  // Calculate gradient of pressure
  mesh->cub_grad_with_central_flux(pr, dpdx, dpdy);

  if(false) {
    // Calculate new velocity intermediate values
    op_par_loop(ins_pressure_update_2d, "ins_pressure_update_2d", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));
  } else {
    // Calculate new velocity intermediate values
    op_par_loop(project_2d_0, "project_2d_0", mesh->cells,
                op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->J, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(proj_rhs_x, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(proj_rhs_y, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

    mesh->mass(proj_rhs_x);
    mesh->mass(proj_rhs_y);

    DG_FP factor = dt * 1.0;
    // DG_FP factor = dt / Cr;
    // op_printf("Cr: %g\n", Cr);
    op_par_loop(project_2d_pen, "project_2d_pen", mesh->cells,
                op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_h, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen, -1, OP_ID, 1, DG_FP_STR, OP_WRITE));

    // Do the 2 vector linear solve with project_mat as the matrix
    int num_cells = 0;
    int num_converge = 0;
    DG_FP num_iter = 0.0;
    op_par_loop(project_2d_cg, "project_2d_cg", mesh->cells,
                op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->J, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_op_xx, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_op_yy, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_op_yx, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_op_xy, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_pen, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(proj_rhs_x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(proj_rhs_y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
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
  }
  */
}

bool INSSolver2D::viscosity() {
  /*
  timer->startTimer("INSSolver2D - Viscosity RHS");
  DG_FP time_n1 = time + dt;

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(visBC[0], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visBC[1], -1, OP_ID, DG_G_NP, DG_FP_STR, OP_WRITE));

  // Get BCs for viscosity solve
  if(mesh->bface2cells) {
    op_par_loop(ins_vis_bc_2d, "ins_vis_bc_2d", mesh->bfaces,
                op_arg_gbl(&time_n1, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->gauss->x,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->y,  0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->nx, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->gauss->ny, 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_READ),
                op_arg_dat(visBC[0], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC),
                op_arg_dat(visBC[1], 0, mesh->bface2cells, DG_G_NP, DG_FP_STR, OP_INC));
  }
  // Set up RHS for viscosity solve
  op_par_loop(ins_vis_copy_2d, "ins_vis_copy_2d", mesh->cells,
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visRHS[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

   mesh->mass(visRHS[0]);
   mesh->mass(visRHS[1]);

  DG_FP factor = reynolds / dt;
  op_par_loop(ins_vis_rhs_2d, "ins_vis_rhs_2d", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(visRHS[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  timer->endTimer("INSSolver2D - Viscosity RHS");

  // Call PETSc linear solver
  timer->startTimer("INSSolver2D - Viscosity Linear Solve");
  factor = g0 * reynolds / dt;
  if(factor != viscosityMatrix->get_factor()) {
    viscosityMatrix->set_factor(factor);
    viscosityMatrix->set_bc_types(vis_bc_types);
    // viscosityMatrix->calc_mat();
    viscositySolver->setFactor(1.0 / factor);
  }
  viscositySolver->set_bcs(visBC[0]);
  bool convergedX = viscositySolver->solve(visRHS[0], vel[(currentInd + 1) % 2][0]);

  viscositySolver->set_bcs(visBC[1]);
  bool convergedY = viscositySolver->solve(visRHS[1], vel[(currentInd + 1) % 2][1]);
  timer->endTimer("INSSolver2D - Viscosity Linear Solve");

  // timer->startTimer("Filtering");
  // filter(mesh, Q[(currentInd + 1) % 2][0]);
  // filter(mesh, Q[(currentInd + 1) % 2][1]);
  // timer->endTimer("Filtering");

  return convergedX && convergedY;
  */
  return true;
}

DG_FP INSSolver2D::get_time() {
  return time;
}

DG_FP INSSolver2D::get_dt() {
  return dt;
}

void INSSolver2D::dump_data(const std::string &filename) {
  timer->startTimer("INSSolver2D - Dump Data");
  op_fetch_data_hdf5_file(mesh->x, filename.c_str());
  op_fetch_data_hdf5_file(mesh->y, filename.c_str());
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
  op_fetch_data_hdf5_file(velT[0], filename.c_str());
  op_fetch_data_hdf5_file(velT[1], filename.c_str());
  op_fetch_data_hdf5_file(velTT[0], filename.c_str());
  op_fetch_data_hdf5_file(velTT[1], filename.c_str());
  op_fetch_data_hdf5_file(pr, filename.c_str());
  op_fetch_data_hdf5_file(mesh->order, filename.c_str());
  timer->endTimer("INSSolver2D - Dump Data");
}
