#include "solvers/3d/ins_solver.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"

extern DGConstants *constants;

#include "timing.h"
#include "config.h"
#include "dg_linear_solvers/petsc_amg.h"
#include "dg_linear_solvers/initial_guess_extrapolation.h"
#include "dg_matrices/3d/factor_mm_poisson_matrix_free_diag_3d.h"
#include "dg_matrices/3d/factor_mm_poisson_matrix_free_block_diag_3d.h"
#include "dg_linear_solvers/petsc_jacobi.h"
#include "dg_linear_solvers/petsc_block_jacobi.h"
#include "dg_dat_pool.h"
#include "dg_utils.h"
#include "dg_abort.h"

#include <string>
#include <iostream>

extern Timing *timer;
extern Config *config;
extern DGDatPool *dg_dat_pool;

INSSolver3D::INSSolver3D(DGMesh3D *m, const DG_FP re) : INSSolverBase3D(m) {
  resuming = false;

  setup_common();

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;

  dt = 0.0;

  reynolds = re;
}

INSSolver3D::INSSolver3D(DGMesh3D *m, const DG_FP re, const std::string &filename, const int iter) : INSSolverBase3D(m, filename) {
  resuming = true;

  setup_common();

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

  reynolds = re;
}

void INSSolver3D::setup_common() {
  tmp_npf_bc = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_tmp_npf_bc");

  tmp_bc_1 = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_tmp_bc_1");

  pr_bc_types  = tmp_bc_1;
  vis_bc_types = tmp_bc_1;
  pr_bc  = tmp_npf_bc;
  vis_bc = tmp_npf_bc;

  // Pressure matrix and solver
  std::string pr_solver = "p-multigrid";
  config->getStr("pressure-solve", "preconditioner", pr_solver);
  pressureSolverType = set_solver_type(pr_solver);
  if(pressureSolverType != LinearSolver::PETSC_PMULTIGRID)
    dg_abort("Only \'p-multigrid\' preconditioner is supported for 3D single phase flow.");
  int pr_over_int = 0;
  config->getInt("pressure-solve", "over_int", pr_over_int);
  if(pr_over_int != 0)
    dg_abort("Cannot over integrate the pressure solve for 3D single phase flow");
  pressureCoarseMatrix = new PoissonCoarseMatrix3D(mesh);
  pressureMatrix = new PoissonMatrixFreeDiag3D(mesh);
  PETScPMultigrid *tmp_pressureSolver = new PETScPMultigrid(mesh);
  tmp_pressureSolver->set_coarse_matrix(pressureCoarseMatrix);
  pressureSolver = tmp_pressureSolver;

  // Viscous matrix and solver
  std::string vis_solver = "inv-mass";
  config->getStr("viscous-solve", "preconditioner", vis_solver);
  viscositySolverType = set_solver_type(vis_solver);
  if(shock_capturing) {
    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      viscosityMatrix = new FactorMMPoissonMatrixFreeDiag3D(mesh);
      viscositySolver = new PETScJacobiSolver(mesh);
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      viscosityMatrix = new FactorMMPoissonMatrixFreeBlockDiag3D(mesh);
      viscositySolver = new PETScBlockJacobiSolver(mesh);
    } else {
      dg_abort("Only \'jacobi\' and \'block-jacobi\' preconditioner is supported for 3D single phase flow with shock capturing.");
    }
  } else {
    if(viscositySolverType != LinearSolver::PETSC_INV_MASS)
      dg_abort("Only \'inv-mass\' preconditioner is supported for 3D single phase flow without shock capturing.");
    viscosityMatrix = new MMPoissonMatrixFree3D(mesh);
    viscositySolver = new PETScInvMassSolver(mesh);
  }

  pressureSolver->set_matrix(pressureMatrix);
  viscositySolver->set_matrix(viscosityMatrix);

  setup_pressure_viscous_solvers(pressureSolver, viscositySolver);
}

INSSolver3D::~INSSolver3D() {
  delete pressureCoarseMatrix;
  delete pressureMatrix;
  delete viscosityMatrix;
  delete pressureSolver;
  delete viscositySolver;
}

void INSSolver3D::init() {
  timer->startTimer("INSSolver3D - Init");
  INSSolverBase3D::init();

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

  double tmp_dt = -1.0;
  config->getDouble("solver-options", "dt", tmp_dt);
  forced_dt = tmp_dt > 0.0;
  if(forced_dt) {
    dt = tmp_dt;
  } else {
    sub_cycle_dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1) * max_vel());
    dt = sub_cycle_dt;
    if(resuming && it_pre_sub_cycle <= 0 && sub_cycles > 1)
      dt = sub_cycle_dt * sub_cycles;
  }
  // dt *= 1e-2;
  op_printf("INS dt is %g\n", dt);

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_bc_types, "ins_3d_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -3, mesh->bface2nodes, 3, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }

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

  timer->endTimer("INSSolver3D - Init");
}

void INSSolver3D::step() {
  timer->startTimer("INSSolver3D - Advection");
  advection();
  timer->endTimer("INSSolver3D - Advection");

  timer->startTimer("INSSolver3D - Filtering");
  if(filter_advec) {
    filter(velT[0]);
    filter(velT[1]);
    filter(velT[2]);
  }
  timer->endTimer("INSSolver3D - Filtering");

  timer->startTimer("INSSolver3D - Pressure");
  pressure();
  timer->endTimer("INSSolver3D - Pressure");

  timer->startTimer("INSSolver3D - Viscosity");
  viscosity();
  timer->endTimer("INSSolver3D - Viscosity");

  currentInd = (currentInd + 1) % 2;
  update_time();
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
}

void INSSolver3D::advection() {
  if(time == 0.0 || sub_cycles < 1 || it_pre_sub_cycle != 0) {
    advec_standard();
  } else {
    advec_sub_cycle();
  }
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
    op_par_loop(ins_3d_pr_dirichlet, "ins_3d_pr_dirichlet", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
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
    dg_abort("\nPressure solve failed to converge\n");

  if(extrapolate_initial_guess)
    add_to_pr_history();
  timer->endTimer("INSSolver3D - Pressure Linear Solve");

  dg_dat_pool->releaseTempDatCells(divVelT);

  zero_dat(dPdN[(currentInd + 1) % 2]);

  timer->startTimer("INSSolver3D - Pressure Projection");
  DGTempDat dpdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdy = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdz = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->grad_with_central_flux(pr, dpdx.dat, dpdy.dat, dpdz.dat);

  project_velocity(dpdx.dat, dpdy.dat, dpdz.dat);

  dg_dat_pool->releaseTempDatCells(dpdx);
  dg_dat_pool->releaseTempDatCells(dpdy);
  dg_dat_pool->releaseTempDatCells(dpdz);
  timer->endTimer("INSSolver3D - Pressure Projection");
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
  DGTempDat tmp_art_vis, tmp_mm_factor;
  viscosityMatrix->set_bc_types(vis_bc_types);
  if(shock_capturing) {
    tmp_art_vis = dg_dat_pool->requestTempDatCells(DG_NP);
    tmp_mm_factor = dg_dat_pool->requestTempDatCells(DG_NP);
    calc_art_vis(velTT[0], tmp_art_vis.dat);

    op_par_loop(add_one, "add_one", mesh->cells,
                op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
    
    op_par_loop(set_val_dg_np, "set_val_dg_np", mesh->cells,
                op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_mm_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
      tmpMatrix->set_mm_factor(tmp_mm_factor.dat);
      tmpMatrix->calc_mat_partial();
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
      tmpMatrix->set_mm_factor(tmp_mm_factor.dat);
      tmpMatrix->calc_mat_partial();
    }
  } else {
    MMPoissonMatrixFree3D *tmpMatrix = dynamic_cast<MMPoissonMatrixFree3D*>(viscosityMatrix);
    PETScInvMassSolver *tmpSolver = dynamic_cast<PETScInvMassSolver*>(viscositySolver);
    if(factor != tmpMatrix->get_factor()) {
      tmpMatrix->set_factor(factor);
      tmpSolver->setFactor(1.0 / factor);
    }
  }

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

  if(shock_capturing) {
    calc_art_vis(velTT[1], tmp_art_vis.dat);

    op_par_loop(add_one, "add_one", mesh->cells,
                op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
      tmpMatrix->calc_mat_partial();
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
      tmpMatrix->calc_mat_partial();
    }
  }

  viscositySolver->set_bcs(vis_bc);
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

  if(shock_capturing) {
    calc_art_vis(velTT[2], tmp_art_vis.dat);

    op_par_loop(add_one, "add_one", mesh->cells,
                op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
      tmpMatrix->calc_mat_partial();
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
      tmpMatrix->calc_mat_partial();
    }
  }

  viscositySolver->set_bcs(vis_bc);
  bool convergedZ = viscositySolver->solve(visRHS[2].dat, vel[(currentInd + 1) % 2][2]);
  if(!convergedZ)
    dg_abort("\nViscosity Z solve failed to converge\n");

  dg_dat_pool->releaseTempDatCells(visRHS[0]);
  dg_dat_pool->releaseTempDatCells(visRHS[1]);
  dg_dat_pool->releaseTempDatCells(visRHS[2]);

  if(shock_capturing) {
    dg_dat_pool->releaseTempDatCells(tmp_art_vis);
    dg_dat_pool->releaseTempDatCells(tmp_mm_factor);
  }
  timer->endTimer("INSSolver3D - Viscosity Linear Solve");
}
