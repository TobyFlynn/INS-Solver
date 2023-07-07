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

INSSolver3D::INSSolver3D(DGMesh3D *m) : INSSolverBase3D(m) {
  resuming = false;

  setup_common();

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;

  dt = 0.0;
  time = 0.0;

  currentInd = 0;
}

INSSolver3D::INSSolver3D(DGMesh3D *m, const std::string &filename, const int iter) : INSSolverBase3D(m, filename) {
  resuming = true;

  setup_common();

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
  tmp_npf_bc = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_tmp_npf_bc");

  tmp_bc_1 = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_tmp_bc_1");

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

  INSSolverBase3D::init(re, refVel);

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

  if(!resuming) {
    op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
                op_arg_dat(dPdN[0], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));

    op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
                op_arg_dat(dPdN[1], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));

    op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                op_arg_dat(pr, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    op_par_loop(zero_np_3, "zero_np_3", mesh->cells,
                op_arg_dat(n[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(n[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(n[0][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    op_par_loop(zero_np_3, "zero_np_3", mesh->cells,
                op_arg_dat(n[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(n[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(n[1][2], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
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

  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, 4 * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("INSSolver3D - Pressure Projection");
  project_velocity();
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
