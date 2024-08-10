#include "solvers/2d/ins_solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "dg_abort.h"

#include "timing.h"
#include "config.h"

#include "dg_linear_solvers/petsc_pmultigrid.h"
#include "dg_linear_solvers/petsc_jacobi.h"
#include "dg_linear_solvers/petsc_block_jacobi.h"
#include "dg_matrices/2d/factor_mm_poisson_matrix_free_diag_2d.h"
#include "dg_matrices/2d/factor_mm_poisson_matrix_free_block_diag_2d.h"

#include "slip_matrix/2d/viscous_matrix.h"
#include "slip_matrix/2d/factor_viscous_matrix.h"

extern Timing *timer;
extern Config *config;
extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

using namespace std;

INSSolver2D::INSSolver2D(DGMesh2D *m, const DG_FP re) : INSSolverBase2D(m) {
  resuming = false;

  setup_common();

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;

  reynolds = re;
}

INSSolver2D::INSSolver2D(DGMesh2D *m, const DG_FP re, const std::string &filename, const int iter) : INSSolverBase2D(m, filename) {
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

void INSSolver2D::setup_common() {
  int tmp_vis = 1;
  config->getInt("solver-options", "viscosity", tmp_vis);
  vis_solve = tmp_vis != 0;
  double tmp_dt = -1.0;
  config->getDouble("solver-options", "force_dt", tmp_dt);
  dt_forced = tmp_dt > 0.0;
  if(dt_forced) dt = tmp_dt;
  int tmp_slip_bcs = 0;
  config->getInt("solver-options", "uses_slip_bcs", tmp_slip_bcs);
  uses_slip_bcs = tmp_slip_bcs != 0;

  // Pressure matrix and solver
  std::string pr_solver = "p-multigrid";
  config->getStr("pressure-solve", "preconditioner", pr_solver);
  pressureSolverType = set_solver_type(pr_solver);
  if(pressureSolverType != LinearSolver::PETSC_PMULTIGRID)
    dg_abort("Only \'p-multigrid\' preconditioner is supported for 2D single phase flow.");
  int pr_over_int = 0;
  config->getInt("pressure-solve", "over_int", pr_over_int);
  if(pr_over_int != 0)
    dg_abort("Cannot over integrate the pressure solve for 2D single phase flow");
  pressureMatrix = new PoissonMatrixFreeDiag2D(mesh);
  pressureCoarseMatrix = new PoissonCoarseMatrix2D(mesh);
  PETScPMultigrid *tmp_pressureSolver = new PETScPMultigrid(mesh);
  tmp_pressureSolver->set_coarse_matrix(pressureCoarseMatrix);
  pressureSolver = tmp_pressureSolver;
  pressureSolver->set_matrix(pressureMatrix);

  // Viscous matrix and solver
  if(uses_slip_bcs) {
    slipViscousSolver = new ViscousSolver2D(mesh);
    if(shock_capturing) {
      slipViscousMatrix = new FactorViscousMatrix2D(mesh);
      slipViscousSolver->set_preconditioner(ViscousSolver2D::FACTOR_INV_MASS);
    } else {
      slipViscousMatrix = new ViscousMatrix2D(mesh);
      slipViscousSolver->set_preconditioner(ViscousSolver2D::FACTOR_INV_MASS);
    }
    slipViscousSolver->set_matrix(slipViscousMatrix);
    double vis_rtol = 1e-8;
    double vis_atol = 1e-9;
    int vis_max_iter = 5000;
    config->getDouble("viscous-solve", "r_tol", vis_rtol);
    config->getDouble("viscous-solve", "a_tol", vis_atol);
    config->getInt("viscous-solve", "max_iter", vis_max_iter);
    slipViscousSolver->set_tol_and_iter(vis_rtol, vis_atol, vis_max_iter);
    bc_data_2 = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_bc_data_2");
    vis_bc_types_2 = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_vis_bc_types_2");
  } else {
    std::string vis_solver = "inv-mass";
    config->getStr("viscous-solve", "preconditioner", vis_solver);
    viscositySolverType = set_solver_type(vis_solver);
    if(shock_capturing) {
      if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
        viscosityMatrix = new FactorMMPoissonMatrixFreeDiag2D(mesh);
        viscositySolver = new PETScJacobiSolver(mesh);
      } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
        viscosityMatrix = new FactorMMPoissonMatrixFreeBlockDiag2D(mesh);
        viscositySolver = new PETScBlockJacobiSolver(mesh);
      } else {
        dg_abort("Only \'jacobi\' and \'block-jacobi\' preconditioner is supported for 2D single phase flow with shock capturing.");
      }
    } else {
      if(viscositySolverType != LinearSolver::PETSC_INV_MASS)
        dg_abort("Only \'inv-mass\' preconditioner is supported for 2D single phase flow without shock capturing.");
      viscosityMatrix = new MMPoissonMatrixFree2D(mesh);
      viscositySolver = new PETScInvMassSolver(mesh);
    }
    viscositySolver->set_matrix(viscosityMatrix);
  }

  setup_pressure_viscous_solvers(pressureSolver, viscositySolver);

  pr_bc_types  = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_pr_bc_types");
  vis_bc_types = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_vis_bc_types");
}

INSSolver2D::~INSSolver2D() {
  delete pressureCoarseMatrix;
  delete pressureMatrix;
  delete pressureSolver;
  if(uses_slip_bcs) {
    delete slipViscousMatrix;
    delete slipViscousSolver;
  } else {
    delete viscosityMatrix;
    delete viscositySolver;
  }
}

void INSSolver2D::init() {
  timer->startTimer("INSSolver2D - Init");
  INSSolverBase2D::init();

  // Set initial conditions
  if(!resuming) {
    op_par_loop(ins_2d_set_ic, "ins_2d_set_ic", mesh->cells,
                op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
                op_arg_dat(vel[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

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
    if(resuming && it_pre_sub_cycle <= 0 && sub_cycles > 1)
      dt = sub_cycle_dt * sub_cycles;
  }
  op_printf("dt: %g\n", dt);

  if(mesh->bface2nodes) {
    op_par_loop(ins_bc_types, "ins_bc_types", mesh->bfaces,
                op_arg_dat(mesh->node_coords, -2, mesh->bface2nodes, 2, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_WRITE));

    op_par_loop(ins_2d_set_pr_bc_type, "ins_2d_set_pr_bc_type", mesh->bfaces,
                op_arg_dat(bc_types,    -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(pr_bc_types, -1, OP_ID, 1, "int", OP_WRITE));

    op_par_loop(ins_2d_set_vis_bc_type, "ins_2d_set_vis_bc_type", mesh->bfaces,
                op_arg_dat(bc_types,     -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }

  pressureCoarseMatrix->set_bc_types(pr_bc_types);
  pressureMatrix->set_bc_types(pr_bc_types);
  pressureSolver->init();
  if(!uses_slip_bcs) {
    viscosityMatrix->set_bc_types(vis_bc_types);
    viscositySolver->init();
  }

  timer->endTimer("INSSolver2D - Init");
}

void INSSolver2D::step() {
  timer->startTimer("INSSolver2D - Advection");
  advection();
  timer->endTimer("INSSolver2D - Advection");

  timer->startTimer("INSSolver2D - Pressure");
  pressure();
  timer->endTimer("INSSolver2D - Pressure");

  timer->startTimer("INSSolver2D - Viscosity");
  if(vis_solve) {
    viscosity();
  } else {
    no_viscosity();
  }
  timer->endTimer("INSSolver2D - Viscosity");

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
void INSSolver2D::advection() {
  if(time == 0.0 || sub_cycles < 1 || it_pre_sub_cycle != 0) {
    advec_standard();
  } else {
    advec_sub_cycle();
  }
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
  // This assumes that the boundaries will always be order DG_ORDER
  op_par_loop(ins_pressure_rhs_2d, "ins_pressure_rhs_2d", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(divVelT.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT.dat);
  mesh->mass(divVelT.dat);
  timer->endTimer("INSSolver2D - Pressure RHS");

  // Call PETSc linear solver
  timer->startTimer("INSSolver2D - Pressure Linear Solve");
  pressureSolver->set_bcs(bc_data);
  bool converged = pressureSolver->solve(divVelT.dat, pr);
  if(!converged)
    dg_abort("Pressure solve did not converge");
  dg_dat_pool->releaseTempDatCells(divVelT);
  timer->endTimer("INSSolver2D - Pressure Linear Solve");

  zero_dat(dPdN[(currentInd + 1) % 2]);

  timer->startTimer("INSSolver2D - Pressure Projection");
  DGTempDat dpdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdy = dg_dat_pool->requestTempDatCells(DG_NP);
  // Calculate gradient of pressure
  mesh->grad_with_central_flux(pr, dpdx.dat, dpdy.dat);

  project_velocity(dpdx.dat, dpdy.dat);

  dg_dat_pool->releaseTempDatCells(dpdx);
  dg_dat_pool->releaseTempDatCells(dpdy);
  timer->endTimer("INSSolver2D - Pressure Projection");

  return converged;
}

bool INSSolver2D::viscosity() {
  timer->startTimer("INSSolver2D - Viscosity RHS");
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

  timer->endTimer("INSSolver2D - Viscosity RHS");
  factor = g0 * reynolds / dt;

  if(uses_slip_bcs) {
    // BCs
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
                  op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                  op_arg_dat(bc_data, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));

      op_par_loop(ins_2d_vis_bc_y, "ins_2d_vis_bc_y", mesh->bfaces,
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
                  op_arg_dat(vis_bc_types_2, -1, OP_ID, 1, "int", OP_WRITE),
                  op_arg_dat(bc_data_2, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
    }

    DGTempDat tmp_art_vis, tmp_mm_factor;
    if(shock_capturing) {
      tmp_art_vis = dg_dat_pool->requestTempDatCells(DG_NP);
      tmp_mm_factor = dg_dat_pool->requestTempDatCells(DG_NP);
      calc_art_vis(velTT[0], tmp_art_vis.dat);
      calc_art_vis(velTT[1], tmp_mm_factor.dat);

      op_par_loop(art_vis_2d_max, "art_vis_2d_max", mesh->cells,
                  op_arg_dat(tmp_mm_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

      op_par_loop(add_one, "add_one", mesh->cells,
                  op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
      
      op_par_loop(set_val_dg_np, "set_val_dg_np", mesh->cells,
                  op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_mm_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

      FactorViscousMatrix2D *tmpMatrix = dynamic_cast<FactorViscousMatrix2D*>(slipViscousMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
      tmpMatrix->set_mm_factor(tmp_mm_factor.dat);

      slipViscousSolver->set_inv_mass_factor(1.0 / factor);
    } else {
      ViscousMatrix2D *tmpMatrix = dynamic_cast<ViscousMatrix2D*>(slipViscousMatrix);
      tmpMatrix->set_factor(factor);

      slipViscousSolver->set_inv_mass_factor(1.0 / factor);
    }

    slipViscousMatrix->set_bc_types(vis_bc_types, vis_bc_types_2);
    slipViscousSolver->set_bcs(bc_data, bc_data_2);

    bool converged = slipViscousSolver->solve(visRHS[0].dat, visRHS[1].dat, vel[(currentInd + 1) % 2][0], vel[(currentInd + 1) % 2][1]);

    if(!converged)
      dg_abort("Viscosity solve did not converge");

    dg_dat_pool->releaseTempDatCells(visRHS[0]);
    dg_dat_pool->releaseTempDatCells(visRHS[1]);
    if(shock_capturing) {
      dg_dat_pool->releaseTempDatCells(tmp_art_vis);
      dg_dat_pool->releaseTempDatCells(tmp_mm_factor);
    }

    return converged;
  } else {
    // Call PETSc linear solver
    timer->startTimer("INSSolver2D - Viscosity Linear Solve");
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
                  op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                  op_arg_dat(bc_data, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
    }
    viscosityMatrix->set_bc_types(vis_bc_types);
    viscositySolver->set_bcs(bc_data);

    DGTempDat tmp_art_vis, tmp_mm_factor;
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
        FactorMMPoissonMatrixFreeDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(tmp_art_vis.dat);
        tmpMatrix->set_mm_factor(tmp_mm_factor.dat);
        tmpMatrix->calc_mat_partial();
      } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
        FactorMMPoissonMatrixFreeBlockDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(tmp_art_vis.dat);
        tmpMatrix->set_mm_factor(tmp_mm_factor.dat);
        tmpMatrix->calc_mat_partial();
      }
    } else {
      MMPoissonMatrixFree2D *tmpMatrix = dynamic_cast<MMPoissonMatrixFree2D*>(viscosityMatrix);
      PETScInvMassSolver *tmpSolver = dynamic_cast<PETScInvMassSolver*>(viscositySolver);
      if(factor != tmpMatrix->get_factor()) {
        tmpMatrix->set_factor(factor);
        tmpSolver->setFactor(1.0 / factor);
      }
    }

    bool convergedX = viscositySolver->solve(visRHS[0].dat, vel[(currentInd + 1) % 2][0]);
    if(!convergedX)
      dg_abort("Viscosity X solve did not converge");

    // Get BCs for viscosity solve
    if(mesh->bface2cells) {
      op_par_loop(ins_2d_vis_bc_y, "ins_2d_vis_bc_y", mesh->bfaces,
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
                  op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE),
                  op_arg_dat(bc_data, -1, OP_ID, DG_NPF, DG_FP_STR, OP_WRITE));
    }
    viscosityMatrix->set_bc_types(vis_bc_types);
    viscositySolver->set_bcs(bc_data);

    if(shock_capturing) {
      calc_art_vis(velTT[1], tmp_art_vis.dat);

      op_par_loop(add_one, "add_one", mesh->cells,
                  op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

      if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
        FactorMMPoissonMatrixFreeDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(tmp_art_vis.dat);
        tmpMatrix->calc_mat_partial();
      } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
        FactorMMPoissonMatrixFreeBlockDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(tmp_art_vis.dat);
        tmpMatrix->calc_mat_partial();
      }
    }

    bool convergedY = viscositySolver->solve(visRHS[1].dat, vel[(currentInd + 1) % 2][1]);
    if(!convergedY)
      dg_abort("Viscosity Y solve did not converge");
    timer->endTimer("INSSolver2D - Viscosity Linear Solve");

    dg_dat_pool->releaseTempDatCells(visRHS[0]);
    dg_dat_pool->releaseTempDatCells(visRHS[1]);

    if(shock_capturing) {
      dg_dat_pool->releaseTempDatCells(tmp_art_vis);
      dg_dat_pool->releaseTempDatCells(tmp_mm_factor);
    }

    return convergedX && convergedY;
  }
}

void INSSolver2D::no_viscosity() {
  op_par_loop(ins_no_vis_2d, "ins_vis_rhs_2d", mesh->cells,
              op_arg_gbl(&g0, 1, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(velTT[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vel[(currentInd + 1) % 2][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(vel[(currentInd + 1) % 2][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
}
