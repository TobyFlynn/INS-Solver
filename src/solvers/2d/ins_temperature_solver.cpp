#include "solvers/2d/ins_temperature_solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"

#include "timing.h"
#include "config.h"
#include "dg_linear_solvers/petsc_amg.h"
#include "dg_linear_solvers/petsc_block_jacobi.h"
#include "dg_linear_solvers/petsc_pmultigrid.h"

extern Timing *timer;
extern Config *config;
extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

using namespace std;

/*********************************************************************************************
 * Temperature Advection Diffusion Solver class that extends the base AdvecDiff Solver class *
 *********************************************************************************************/
TemperatureAdvecDiffSolver2D::TemperatureAdvecDiffSolver2D(DGMesh2D *m) : AdvecDiffSolver2D(m) {}

void TemperatureAdvecDiffSolver2D::set_bc_types(op_dat bc) {
  bc_types = bc;
}

void TemperatureAdvecDiffSolver2D::bc_kernel_advec(op_dat val, op_dat u, op_dat v, op_dat out) {
  op_par_loop(temperature_advec_2d_bc, "temperature_advec_2d_bc", mesh->bfaces,
              op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
}

void TemperatureAdvecDiffSolver2D::bc_kernel_diff_0(op_dat val, op_dat out_x, op_dat out_y) {
  op_par_loop(temperature_diff_2d_bc_0, "temperature_diff_2d_bc_0", mesh->bfaces,
              op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out_x, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(out_y, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
}

void TemperatureAdvecDiffSolver2D::bc_kernel_diff_1(op_dat val, op_dat val_x, op_dat val_y, op_dat vis, op_dat out) {
  op_par_loop(temperature_diff_2d_bc_1, "temperature_diff_2d_bc_1", mesh->bfaces,
              op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val_x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val_y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vis, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
}

/************************
 * Main INS Temperature Solver class *
 ************************/
INSTemperatureSolver2D::INSTemperatureSolver2D(DGMesh2D *m) : INSSolverBase2D(m) {
  resuming = false;

  setup_common();

  currentInd = 0;

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;
}

INSTemperatureSolver2D::INSTemperatureSolver2D(DGMesh2D *m, const std::string &filename, const int iter) : INSSolverBase2D(m, filename) {
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

void INSTemperatureSolver2D::setup_common() {
  double tmp_dt = -1.0;
  config->getDouble("solver-options", "force_dt", tmp_dt);
  dt_forced = tmp_dt > 0.0;
  if(dt_forced) dt = tmp_dt;

  // Pressure matrix and solver
  std::string pr_solver = "p-multigrid";
  config->getStr("pressure-solve", "preconditioner", pr_solver);
  if(pr_solver != "p-multigrid") 
    throw std::runtime_error("Only \'p-multigrid\' preconditioner is supported for 2D temperature + single phase flow.");
  pressureMatrix = new FactorPoissonMatrixFreeDiag2D(mesh);
  pressureCoarseMatrix = new FactorPoissonCoarseMatrix2D(mesh);
  pressureSolver = new PETScPMultigrid(mesh);
  pressureSolver->set_coarse_matrix(pressureCoarseMatrix);

  // Viscous matrix and solver
  std::string vis_solver = "inv-mass";
  config->getStr("viscous-solve", "preconditioner", vis_solver);
  if(vis_solver != "inv-mass") 
    throw std::runtime_error("Only \'inv-mass\' preconditioner is supported for 2D temperature + single phase flow.");
  viscosityMatrix = new MMPoissonMatrixFree2D(mesh);
  viscositySolver = new PETScInvMassSolver(mesh);

  advecDiffSolver = new TemperatureAdvecDiffSolver2D(mesh);

  pressureSolver->set_matrix(pressureMatrix);
  viscositySolver->set_matrix(viscosityMatrix);

  setup_pressure_viscous_solvers(pressureSolver, viscositySolver);

  pr_bc_types  = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_pr_bc_types");
  vis_bc_types = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_vis_bc_types");

  temperature = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_temperature");
  rho = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_rho");
}

INSTemperatureSolver2D::~INSTemperatureSolver2D() {
  delete pressureCoarseMatrix;
  delete pressureMatrix;
  delete viscosityMatrix;
  delete pressureSolver;
  delete viscositySolver;
  delete advecDiffSolver;
}

void INSTemperatureSolver2D::init(const DG_FP re, const DG_FP refVel) {
  timer->startTimer("INSTemperatureSolver2D - Init");
  INSSolverBase2D::init(re, refVel);

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

    op_par_loop(ins_2d_set_ic_temperature, "ins_2d_set_ic", mesh->cells,
                op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(temperature, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    update_rho();

    zero_dat(dPdN[0]);
    zero_dat(dPdN[1]);
    zero_dat(n[0][0]);
    zero_dat(n[0][1]);
    zero_dat(n[1][0]);
    zero_dat(n[1][1]);
    zero_dat(pr);
  }

  sub_cycle_dt = h / (DG_ORDER * DG_ORDER * max_vel());
  if(!dt_forced) {
    dt = sub_cycle_dt;
    if(resuming)
      dt = sub_cycles > 1 ? sub_cycle_dt * sub_cycles : sub_cycle_dt;
  }
  op_printf("dt: %g\n", dt);

  time = dt * currentInd;
  currentInd = currentInd % 2;

  if(mesh->bface2nodes) {
    op_par_loop(ins_bc_types, "ins_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -2, mesh->bface2nodes, 2, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_WRITE));

    op_par_loop(ins_2d_set_pr_bc_type, "ins_2d_set_pr_bc_type", mesh->bfaces,
                op_arg_dat(bc_types,    -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(pr_bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }
  advecDiffSolver->set_bc_types(pr_bc_types);

  pressureCoarseMatrix->set_bc_types(pr_bc_types);
  pressureMatrix->set_bc_types(pr_bc_types);
  pressureSolver->init();
  viscositySolver->init();

  timer->endTimer("INSTemperatureSolver2D - Init");
}

void INSTemperatureSolver2D::step() {
  timer->startTimer("INSTemperatureSolver2D - Advection");
  advection();
  timer->endTimer("INSTemperatureSolver2D - Advection");

  timer->startTimer("INSTemperatureSolver2D - Pressure");
  pressure();
  timer->endTimer("INSTemperatureSolver2D - Pressure");

  timer->startTimer("INSTemperatureSolver2D - Viscosity");
  viscosity();
  timer->endTimer("INSTemperatureSolver2D - Viscosity");

  timer->startTimer("INSTemperatureSolver2D - Update T");
  update_temperature();
  timer->endTimer("INSTemperatureSolver2D - Update T");

  currentInd = (currentInd + 1) % 2;
  update_time();

  g0 = 1.5;
  a0 = 2.0;
  a1 = -0.5;
  b0 = 2.0;
  b1 = -1.0;

  if(!dt_forced) {
    if(it_pre_sub_cycle > 1) {
      it_pre_sub_cycle--;
    } else {
      sub_cycle_dt = h / (DG_ORDER * DG_ORDER * max_vel());
      dt = sub_cycles > 1 ? sub_cycle_dt * sub_cycles : sub_cycle_dt;
      it_pre_sub_cycle = 0;
    }
  }
}

// Calculate Nonlinear Terms
void INSTemperatureSolver2D::advection() {
  if(time == 0.0 || sub_cycles < 1 || it_pre_sub_cycle != 0) {
    advec_standard();
  } else {
    advec_sub_cycle();
  }
}

bool INSTemperatureSolver2D::pressure() {
  timer->startTimer("INSTemperatureSolver2D - Pressure RHS");
  DGTempDat divVelT = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat curlVel = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat gradCurlVel[2];
  gradCurlVel[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  gradCurlVel[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->div_with_central_flux(velT[0], velT[1], divVelT.dat);

  // Calculate RHS of pressure solve
  const DG_FP div_factor = -1.0 / dt;
  op_par_loop(mp_ins_2d_pr_div_factor, "mp_ins_2d_pr_div_factor", mesh->cells,
              op_arg_gbl(&div_factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(divVelT.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  mesh->curl(vel[currentInd][0], vel[currentInd][1], curlVel.dat);
  mesh->grad(curlVel.dat, gradCurlVel[0].dat, gradCurlVel[1].dat);

  dg_dat_pool->releaseTempDatCells(curlVel);

  // Apply Neumann pressure boundary conditions
  if(mesh->bface2cells) {
    op_par_loop(mp_ins_2d_pr_0, "mp_ins_2d_pr_0", mesh->bfaces,
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
                op_arg_dat(rho, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dPdN[currentInd],   0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

  dg_dat_pool->releaseTempDatCells(gradCurlVel[0]);
  dg_dat_pool->releaseTempDatCells(gradCurlVel[1]);

  // Apply Dirichlet BCs
  if(mesh->bface2cells) {
    op_par_loop(ins_2d_pr_bc, "ins_2d_pr_bc", mesh->bfaces,
                op_arg_gbl(&time, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(pr_bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(bc_data, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  // Calculate RHS of pressure solve
  op_par_loop(mp_ins_2d_pr_bc_1, "mp_ins_2d_pr_bc_1", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT.dat);
  mesh->mass(divVelT.dat);

  DGTempDat pr_factor = dg_dat_pool->requestTempDatCells(DG_NP);
  op_par_loop(reciprocal, "reciprocal", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pr_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("INSTemperatureSolver2D - Pressure RHS");

  // Call PETSc linear solver
  timer->startTimer("INSTemperatureSolver2D - Pressure Linear Solve");
  pressureMatrix->set_factor(pr_factor.dat);
  pressureCoarseMatrix->set_factor(pr_factor.dat);
  pressureMatrix->set_bc_types(pr_bc_types);
  pressureCoarseMatrix->set_bc_types(pr_bc_types);
  pressureSolver->set_coarse_matrix(pressureCoarseMatrix);
  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_bcs(bc_data);
  bool converged = pressureSolver->solve(divVelT.dat, pr);
  if(!converged)
    throw std::runtime_error("Pressure solve did not converge");

  dg_dat_pool->releaseTempDatCells(pr_factor);
  dg_dat_pool->releaseTempDatCells(divVelT);
  timer->endTimer("INSTemperatureSolver2D - Pressure Linear Solve");

  zero_dat(dPdN[(currentInd + 1) % 2]);

  timer->startTimer("INSTemperatureSolver2D - Pressure Projection");
  DGTempDat dpdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdy = dg_dat_pool->requestTempDatCells(DG_NP);
  // Calculate gradient of pressure
  mesh->grad_over_int_with_central_flux(pr, dpdx.dat, dpdy.dat);

  op_par_loop(mp_ins_2d_pr_2, "mp_ins_2d_pr_2", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdx.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  project_velocity(dpdx.dat, dpdy.dat);

  shock_capture(velTT[0], velTT[1]);

  dg_dat_pool->releaseTempDatCells(dpdx);
  dg_dat_pool->releaseTempDatCells(dpdy);
  timer->endTimer("INSTemperatureSolver2D - Pressure Projection");

  return converged;
}

bool INSTemperatureSolver2D::viscosity() {
  timer->startTimer("INSTemperatureSolver2D - Viscosity RHS");
  DG_FP time_n1 = time + dt;

  // Set up RHS for viscosity solve
  DGTempDat visRHS[2];
  visRHS[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  visRHS[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  op_par_loop(ins_vis_copy_2d, "ins_vis_copy_2d", mesh->cells,
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(visRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

   mesh->mass(visRHS[0].dat);
   mesh->mass(visRHS[1].dat);

  DG_FP factor = reynolds / dt;
  op_par_loop(ins_vis_rhs_2d, "ins_vis_rhs_2d", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(visRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  timer->endTimer("INSTemperatureSolver2D - Viscosity RHS");

  // Call PETSc linear solver
  timer->startTimer("INSTemperatureSolver2D - Viscosity Linear Solve");
  op_par_loop(ins_2d_set_vis_x_bc_type, "ins_2d_set_vis_x_bc_type", mesh->bfaces,
              op_arg_dat(bc_types,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  viscosityMatrix->set_bc_types(vis_bc_types);

  factor = g0 * reynolds / dt;
  if(factor != viscosityMatrix->get_factor()) {
    viscosityMatrix->set_factor(factor);
    // viscosityMatrix->calc_mat();
    viscositySolver->setFactor(1.0 / factor);
  }

  // Get BCs for viscosity solve
  if(mesh->bface2cells) {
    op_par_loop(ins_2d_vis_bc_x, "ins_2d_vis_bc_x", mesh->bfaces,
                op_arg_gbl(&time_n1, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&g0, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(bc_data, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  viscositySolver->set_bcs(bc_data);
  bool convergedX = viscositySolver->solve(visRHS[0].dat, vel[(currentInd + 1) % 2][0]);
  if(!convergedX)
    throw std::runtime_error("Viscosity X solve did not converge");

  op_par_loop(ins_2d_set_vis_y_bc_type, "ins_2d_set_vis_y_bc_type", mesh->bfaces,
              op_arg_dat(bc_types,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  viscosityMatrix->set_bc_types(vis_bc_types);

  // Get BCs for viscosity solve
  if(mesh->bface2cells) {
    op_par_loop(ins_2d_vis_bc_y, "ins_2d_vis_bc_x", mesh->bfaces,
                op_arg_gbl(&time_n1, 1, DG_FP_STR, OP_READ),
                op_arg_gbl(&g0, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[0], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(velTT[1], 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(bc_data, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
  }

  viscositySolver->set_bcs(bc_data);
  bool convergedY = viscositySolver->solve(visRHS[1].dat, vel[(currentInd + 1) % 2][1]);
  if(!convergedY)
    throw std::runtime_error("Viscosity Y solve did not converge");
  timer->endTimer("INSTemperatureSolver2D - Viscosity Linear Solve");

  dg_dat_pool->releaseTempDatCells(visRHS[0]);
  dg_dat_pool->releaseTempDatCells(visRHS[1]);

  // timer->startTimer("Filtering");
  // filter(mesh, Q[(currentInd + 1) % 2][0]);
  // filter(mesh, Q[(currentInd + 1) % 2][1]);
  // timer->endTimer("Filtering");

  return convergedX && convergedY;
}

void INSTemperatureSolver2D::update_rho() {
  op_par_loop(ins_2d_update_rho_temperature, "ins_2d_update_rho_temperature", mesh->cells,
              op_arg_dat(temperature, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}

void INSTemperatureSolver2D::update_temperature() {
  // Set thermal conductivity / heat capacity per unit volume
  DGTempDat temperature_diff_coef = dg_dat_pool->requestTempDatCells(DG_NP);
  op_par_loop(ins_2d_set_temperature_diff, "ins_2d_set_temperature_diff", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ)
              op_arg_dat(temperature_diff_coef.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  advecDiffSolver->step(temperature, vel[(currentInd + 1) % 2][0], vel[(currentInd + 1) % 2][1], temperature_diff_coef.dat, dt);

  dg_dat_pool->releaseTempDatCells(temperature_diff_coef);

  update_rho();
}

void INSTemperatureSolver2D::dump_checkpoint_data(const std::string &filename) {
  INSSolverBase2D::dump_checkpoint_data(filename);
  op_fetch_data_hdf5_file(temperature, filename.c_str());

  // TODO save constants in same HDF5 file
}

void INSTemperatureSolver2D::dump_visualisation_data(const std::string &filename) {
  INSSolverBase2D::dump_visualisation_data(filename);

  if(values_to_save.count("temperature") != 0) {
    op_fetch_data_hdf5_file(temperature, filename.c_str());
  }

  if(values_to_save.count("rho") != 0) {
    op_fetch_data_hdf5_file(rho, filename.c_str());
  }
}
