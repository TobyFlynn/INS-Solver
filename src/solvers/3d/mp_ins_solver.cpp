#include "solvers/3d/mp_ins_solver.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_utils.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"
#include "config.h"
#include "op2_utils.h"
#include "dg_linear_solvers/petsc_amg.h"
#include "dg_linear_solvers/petsc_block_jacobi.h"
#include "dg_linear_solvers/initial_guess_extrapolation.h"
#include "dg_abort.h"

#include <string>

extern Timing *timer;
extern DGConstants *constants;
extern Config *config;
extern DGDatPool *dg_dat_pool;

MPINSSolver3D::MPINSSolver3D(DGMesh3D *m) : INSSolverBase3D(m) {
  resuming = false;

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;

  dt = 0.0;
  time = 0.0;

  currentInd = 0;

  lsSolver = new LevelSetSolver3D(mesh);

  setup_common();
}

MPINSSolver3D::MPINSSolver3D(DGMesh3D *m, const std::string &filename, const int iter) : INSSolverBase3D(m, filename) {
  resuming = true;

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

  lsSolver = new LevelSetSolver3D(mesh, filename);

  setup_common();
}

void MPINSSolver3D::setup_common() {
  // Pressure matrix and solver
  std::string pr_solver = "p-multigrid";
  config->getStr("pressure-solve", "preconditioner", pr_solver);
  if(pr_solver != "p-multigrid") 
    throw std::runtime_error("Only \'p-multigrid\' preconditioner is supported for 3D multiphase flow.");
  int pr_over_int = 0;
  config->getInt("pressure-solve", "over_int", pr_over_int);
  if(pr_over_int != 0)
    throw std::runtime_error("Over integrating the pressure solve for 3D multiphase flow is not yet supported");
  coarsePressureMatrix = new FactorPoissonCoarseMatrix3D(mesh);
  pressureMatrix = new FactorPoissonMatrixFreeDiag3D(mesh);
  pressureSolver = new PETScPMultigrid(mesh);

  // Viscous matrix and solver
  std::string vis_solver = "jacobi";
  config->getStr("viscous-solve", "preconditioner", vis_solver);
  if(vis_solver != "jacobi") 
    throw std::runtime_error("Only \'jacobi\' preconditioner is supported for 3D multiphase flow.");
  viscosityMatrix = new FactorMMPoissonMatrixFreeDiag3D(mesh);
  viscositySolver = new PETScJacobiSolver(mesh);
  
  pressureSolver->set_matrix(pressureMatrix);
  viscositySolver->set_matrix(viscosityMatrix);

  setup_pressure_viscous_solvers(pressureSolver, viscositySolver);

  rho = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_rho");
  mu  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_mu");

  tmp_npf_bc = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_tmp_npf_bc");
  tmp_bc_1 = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_tmp_bc_1");
  art_vis  = op_decl_dat(mesh->cells, 1, DG_FP_STR, (DG_FP *)NULL, "ins_solver_art_vis");

  pr_bc_types  = tmp_bc_1;
  vis_bc_types = tmp_bc_1;
  pr_bc  = tmp_npf_bc;
  vis_bc = tmp_npf_bc;
}

MPINSSolver3D::~MPINSSolver3D() {
  delete coarsePressureMatrix;
  delete pressureMatrix;
  delete viscosityMatrix;
  delete pressureSolver;
  delete viscositySolver;
  delete lsSolver;
}

void MPINSSolver3D::init(const DG_FP re, const DG_FP refVel) {
  timer->startTimer("MPINSSolver3D - Init");

  INSSolverBase3D::init(re, refVel);

  reynolds = re;

  lsSolver->init();

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

  sub_cycle_dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1) * max_vel());
  dt = sub_cycle_dt;
  if(resuming)
    dt = sub_cycles > 1 ? sub_cycle_dt * sub_cycles : sub_cycle_dt;
  // dt *= 1e-2;
  op_printf("INS dt is %g\n", dt);
  time = dt * currentInd;
  currentInd = currentInd % 2;

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_bc_types, "ins_3d_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -3, mesh->bface2nodes, 3, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }

  lsSolver->getRhoMu(rho, mu);

  if(!resuming) {
    zero_dat(dPdN[0]);
    zero_dat(dPdN[1]);
    zero_dat(pr);
    zero_dat(n[0][0]);
    zero_dat(n[0][1]);
    zero_dat(n[0][2]);
    zero_dat(n[1][0]);
    zero_dat(n[1][1]);
    zero_dat(n[1][2]);
  }

  pressureSolver->init();
  viscositySolver->init();

  timer->endTimer("MPINSSolver3D - Init");
}

void MPINSSolver3D::step() {
  timer->startTimer("MPINSSolver3D - Advection");
  advection();
  timer->endTimer("MPINSSolver3D - Advection");

  timer->startTimer("MPINSSolver3D - Filtering");
  if(filter_advec) {
    filter(velT[0]);
    filter(velT[1]);
    filter(velT[2]);
  }
  timer->endTimer("MPINSSolver3D - Filtering");

  timer->startTimer("MPINSSolver3D - Pressure");
  pressure();
  timer->endTimer("MPINSSolver3D - Pressure");

  timer->startTimer("MPINSSolver3D - Viscosity");
  viscosity();
  timer->endTimer("MPINSSolver3D - Viscosity");

  timer->startTimer("MPINSSolver3D - Surface");
  surface();
  timer->endTimer("MPINSSolver3D - Surface");

  currentInd = (currentInd + 1) % 2;
  update_time();
  g0 = 1.5;
  a0 = 2.0;
  a1 = -0.5;
  b0 = 2.0;
  b1 = -1.0;
  sub_cycle_dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1) * max_vel());
  if(it_pre_sub_cycle > 1) {
    it_pre_sub_cycle--;
    dt = sub_cycle_dt;
  } else {
    dt = sub_cycles > 1 ? sub_cycle_dt * sub_cycles : sub_cycle_dt;
    it_pre_sub_cycle = 0;
  }
}

void MPINSSolver3D::advection() {
  if(time == 0.0 || sub_cycles < 1 || it_pre_sub_cycle != 0) {
    advec_standard();
  } else {
    advec_sub_cycle();
  }
}

void MPINSSolver3D::pressure() {
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

  op_par_loop(mp_ins_3d_pr_mu, "mp_ins_3d_pr_mu", mesh->cells,
              op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(curlVel[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(curlVel[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(curlVel[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  mesh->curl(curlVel[0].dat, curlVel[1].dat, curlVel[2].dat, curl2Vel[0].dat,
             curl2Vel[1].dat, curl2Vel[2].dat);

  dg_dat_pool->releaseTempDatCells(curlVel[0]);
  dg_dat_pool->releaseTempDatCells(curlVel[1]);
  dg_dat_pool->releaseTempDatCells(curlVel[2]);

  if(mesh->bface2cells) {
    op_par_loop(mp_ins_3d_pr_0, "mp_ins_3d_pr_0", mesh->bfaces,
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
                op_arg_dat(rho, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dPdN[currentInd], 0, mesh->bface2cells, 4 * DG_NPF, DG_FP_STR, OP_INC));
  }

  dg_dat_pool->releaseTempDatCells(curl2Vel[0]);
  dg_dat_pool->releaseTempDatCells(curl2Vel[1]);
  dg_dat_pool->releaseTempDatCells(curl2Vel[2]);

  DGTempDat pr_factor = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(mp_ins_3d_pr_1, "mp_ins_3d_pr_1", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(divVelT.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pr_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT.dat);
  mesh->mass(divVelT.dat);

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_pr_2, "ins_3d_pr_2", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(pr_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(pr_bc, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  timer->startTimer("MPINSSolver3D - Pressure Linear Solve");
  pressureMatrix->set_factor(pr_factor.dat);
  coarsePressureMatrix->set_factor(pr_factor.dat);
  pressureMatrix->set_bc_types(pr_bc_types);
  coarsePressureMatrix->set_bc_types(pr_bc_types);
  pressureSolver->set_coarse_matrix(coarsePressureMatrix);
  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_bcs(pr_bc);

  if(extrapolate_initial_guess)
    initial_guess_extrapolation(mesh, pr_history, pr, time + dt);

  bool converged = pressureSolver->solve(divVelT.dat, pr);
  if(!converged)
    dg_abort("\nPressure solve failed to converge\n");

  if(extrapolate_initial_guess)
    add_to_pr_history();
  timer->endTimer("MPINSSolver3D - Pressure Linear Solve");

  dg_dat_pool->releaseTempDatCells(pr_factor);
  dg_dat_pool->releaseTempDatCells(divVelT);

  zero_dat(dPdN[(currentInd + 1) % 2]);

  timer->startTimer("MPINSSolver3D - Pressure Projection");
  DGTempDat dpdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdy = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdz = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->grad_with_central_flux(pr, dpdx.dat, dpdy.dat, dpdz.dat);

  op_par_loop(mp_ins_3d_pr_2, "mp_ins_3d_pr_2", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdx.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdz.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  project_velocity(dpdx.dat, dpdy.dat, dpdz.dat);

  // shock_capture(velTT[0], velTT[1], velTT[2]);

  dg_dat_pool->releaseTempDatCells(dpdx);
  dg_dat_pool->releaseTempDatCells(dpdy);
  dg_dat_pool->releaseTempDatCells(dpdz);
  timer->endTimer("MPINSSolver3D - Pressure Projection");
}

void MPINSSolver3D::viscosity() {
  DGTempDat visRHS[3];
  visRHS[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  visRHS[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  visRHS[2] = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat vis_mm_factor = dg_dat_pool->requestTempDatCells(DG_NP);

  DG_FP factor  = reynolds / dt;
  DG_FP factor2 = g0 * reynolds / dt;
  op_par_loop(mp_ins_3d_vis_0, "mp_ins_3d_vis_0", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&factor2, 1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mu,  -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(art_vis,  -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visRHS[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vis_mm_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->mass(visRHS[0].dat);
  mesh->mass(visRHS[1].dat);
  mesh->mass(visRHS[2].dat);

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

  timer->startTimer("Vis Linear Solve");
  DGTempDat tmp_art_vis = dg_dat_pool->requestTempDatCells(DG_NP);
  if(shock_capturing) {
    calc_art_vis(velTT[0], tmp_art_vis.dat);

    op_par_loop(mp_ins_3d_add_mu, "mp_ins_3d_add_mu", mesh->cells,
                op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    viscosityMatrix->set_factor(tmp_art_vis.dat);
  } else {
    viscosityMatrix->set_factor(mu);
  }
  viscosityMatrix->set_mm_factor(vis_mm_factor.dat);
  viscosityMatrix->set_bc_types(vis_bc_types);
  viscosityMatrix->calc_mat_partial();
  viscositySolver->set_bcs(vis_bc);
  bool convergedX = viscositySolver->solve(visRHS[0].dat, vel[(currentInd + 1) % 2][0]);
  if(!convergedX)
    dg_abort("\nViscosity X solve failed to converge\n");

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
  if(shock_capturing) {
    calc_art_vis(velTT[1], tmp_art_vis.dat);

    op_par_loop(mp_ins_3d_add_mu, "mp_ins_3d_add_mu", mesh->cells,
                op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    viscosityMatrix->set_factor(tmp_art_vis.dat);
  } else {
    viscosityMatrix->set_factor(mu);
  }
  bool convergedY = viscositySolver->solve(visRHS[1].dat, vel[(currentInd + 1) % 2][1]);
  if(!convergedY)
    dg_abort("\nViscosity Y solve failed to converge\n");

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
  if(shock_capturing) {
    calc_art_vis(velTT[2], tmp_art_vis.dat);

    op_par_loop(mp_ins_3d_add_mu, "mp_ins_3d_add_mu", mesh->cells,
                op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    viscosityMatrix->set_factor(tmp_art_vis.dat);
  } else {
    viscosityMatrix->set_factor(mu);
  }
  bool convergedZ = viscositySolver->solve(visRHS[2].dat, vel[(currentInd + 1) % 2][2]);
  if(!convergedZ)
    dg_abort("\nViscosity Z solve failed to converge\n");

  dg_dat_pool->releaseTempDatCells(tmp_art_vis);
  dg_dat_pool->releaseTempDatCells(vis_mm_factor);
  dg_dat_pool->releaseTempDatCells(visRHS[0]);
  dg_dat_pool->releaseTempDatCells(visRHS[1]);
  dg_dat_pool->releaseTempDatCells(visRHS[2]);
  timer->endTimer("Vis Linear Solve");
}

void MPINSSolver3D::surface() {
  lsSolver->set_bc_types(bc_types);
  const int num_advec_steps = it_pre_sub_cycle != 0 ? 1 : std::max(sub_cycles, 1);
  lsSolver->step(vel[(currentInd + 1) % 2][0], vel[(currentInd + 1) % 2][1], vel[(currentInd + 1) % 2][2],
                 num_advec_steps == 1 ? dt : sub_cycle_dt, num_advec_steps);
  timer->startTimer("MPINSSolver3D - Filtering");
  if(filter_advec) {
    filter(lsSolver->s);
  }
  timer->endTimer("MPINSSolver3D - Filtering");
  lsSolver->getRhoMu(rho, mu);
}

void MPINSSolver3D::dump_checkpoint_data(const std::string &filename) {
  INSSolverBase3D::dump_checkpoint_data(filename);
  op_fetch_data_hdf5_file(lsSolver->s, filename.c_str());

  // TODO save constants in same HDF5 file
}

void MPINSSolver3D::dump_visualisation_data(const std::string &filename) {
  INSSolverBase3D::dump_visualisation_data(filename);

  if(values_to_save.count("mu") != 0) {
    if(shock_capturing) {
      DGTempDat tmp_art_vis = dg_dat_pool->requestTempDatCells(DG_NP);
      calc_art_vis(velTT[0], tmp_art_vis.dat);

      op_par_loop(mp_ins_3d_add_mu, "mp_ins_3d_add_mu", mesh->cells,
                  op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

      op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
                  op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

      dg_dat_pool->releaseTempDatCells(tmp_art_vis);
    }
    op_fetch_data_hdf5_file(mu, filename.c_str());
  }

  if(values_to_save.count("rho") != 0) {
    op_fetch_data_hdf5_file(rho, filename.c_str());
  }

  if(values_to_save.count("level_set") != 0) {
    op_fetch_data_hdf5_file(lsSolver->s, filename.c_str());
  }
}
