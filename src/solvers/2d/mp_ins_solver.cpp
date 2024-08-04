#include "solvers/2d/mp_ins_solver.h"

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

#include "dg_linear_solvers/petsc_jacobi.h"
#include "dg_linear_solvers/petsc_block_jacobi.h"
#include "dg_matrices/2d/factor_mm_poisson_matrix_free_diag_2d.h"
#include "dg_matrices/2d/factor_mm_poisson_matrix_free_block_diag_2d.h"

extern Timing *timer;
extern Config *config;
extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

using namespace std;

MPINSSolver2D::MPINSSolver2D(DGMesh2D *m, const DG_FP re) : INSSolverBase2D(m) {
  resuming = false;

  setup_common();

  lsSolver = new LevelSetSolver2D(m);

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;

  reynolds = re;
}

MPINSSolver2D::MPINSSolver2D(DGMesh2D *m, const DG_FP re, const std::string &filename, const int iter) : INSSolverBase2D(m, filename) {
  resuming = true;

  setup_common();

  lsSolver = new LevelSetSolver2D(mesh, filename);

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

void MPINSSolver2D::setup_common() {
  int tmp_st = 0;
  config->getInt("solver-options", "surface_tension", tmp_st);
  surface_tension = tmp_st == 1;
  int tmp_st_oi = 0;
  config->getInt("solver-options", "over_int_surface_tension", tmp_st_oi);
  over_int_surface_tension = tmp_st_oi == 1;
  double tmp_dt = -1.0;
  config->getDouble("solver-options", "force_dt", tmp_dt);
  dt_forced = tmp_dt > 0.0;
  if(dt_forced) dt = tmp_dt;
  int tmp_grav = 0;
  config->getInt("solver-options", "gravity_modified_pressure", tmp_grav);
  gravity_modified_pressure = tmp_grav == 1;
  int tmp_slip_bcs = 0;
  config->getInt("solver-options", "uses_slip_bcs", tmp_slip_bcs);
  uses_slip_bcs = tmp_slip_bcs != 0;

  if(gravity && gravity_modified_pressure)
    dg_abort("Do not use both \'gravity\' and \'gravity_modified_pressure\'");

  // Pressure matrix and solver
  std::string pr_solver = "p-multigrid";
  config->getStr("pressure-solve", "preconditioner", pr_solver);
  pressureSolverType = set_solver_type(pr_solver);
  if(pressureSolverType != LinearSolver::PETSC_PMULTIGRID)
    dg_abort("Only \'p-multigrid\' preconditioner is supported for 2D multiphase flow.");
  int tmp_pr_over_int = 0;
  config->getInt("pressure-solve", "over_int", tmp_pr_over_int);
  pr_over_int = tmp_pr_over_int != 0;
  if(pr_over_int) {
    pressureMatrix = new FactorPoissonMatrixFreeDiagOI2D(mesh);
  } else {
    pressureMatrix = new FactorPoissonMatrixFreeDiag2D(mesh);
  }
  pressureCoarseMatrix = new FactorPoissonCoarseMatrix2D(mesh);
  pressureSolver = new PETScPMultigrid(mesh);
  pressureSolver->set_coarse_matrix(pressureCoarseMatrix);
  pressureSolver->set_matrix(pressureMatrix);

  // Viscous matrix and solver
  std::string vis_solver = "jacobi";
  config->getStr("viscous-solve", "preconditioner", vis_solver);
  viscositySolverType = set_solver_type(vis_solver);
  if(uses_slip_bcs) {
    slipViscousSolver = new ViscousSolver2D(mesh);
    if(vis_solver == "jacobi") {
      slipViscousMatrix = new FactorViscousMatrix2D(mesh, true, false);
      slipViscousSolver->set_preconditioner(ViscousSolver2D::JACOBI);
    } else if(vis_solver == "block-jacobi") {
      slipViscousMatrix = new FactorViscousMatrix2D(mesh, false, true);
      slipViscousSolver->set_preconditioner(ViscousSolver2D::BLOCK_JACOBI);
    } else if(vis_solver == "inv-mass") {
      slipViscousMatrix = new FactorViscousMatrix2D(mesh, false, false);
      slipViscousSolver->set_preconditioner(ViscousSolver2D::RECP_FACTOR_DAT_INV_MASS);
    } else {
      dg_abort("Unsupported preconditioner for slip BCs");
    }
    slipViscousSolver->set_matrix(slipViscousMatrix);
    slipViscousSolver->set_tol_and_iter(1e-8, 1e-9, 1000);
    bc_data_2 = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_bc_data_2");
    vis_bc_types_2 = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_vis_bc_types_2");
  } else {
    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      viscosityMatrix = new FactorMMPoissonMatrixFreeDiag2D(mesh);
      viscositySolver = new PETScJacobiSolver(mesh);
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      viscosityMatrix = new FactorMMPoissonMatrixFreeBlockDiag2D(mesh);
      viscositySolver = new PETScBlockJacobiSolver(mesh);
    } else {
      dg_abort("Only \'jacobi\' preconditioner is supported for 2D multiphase flow.");
    }
    viscositySolver->set_matrix(viscosityMatrix);
  }

  setup_pressure_viscous_solvers(pressureSolver, viscositySolver);

  rho = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_rho");
  mu  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_mu");
  dPdN_oi[0] = op_decl_dat(mesh->cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_dPdN_oi_0");
  dPdN_oi[1] = op_decl_dat(mesh->cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_dPdN_oi_1");

  if(surface_tension) {
    st[0][0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st00");
    st[0][1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st01");
    st[1][0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st10");
    st[1][1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st11");
  }

  pr_bc_types  = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_pr_bc_types");
  vis_bc_types = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_vis_bc_types");
}

MPINSSolver2D::~MPINSSolver2D() {
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

void MPINSSolver2D::init() {
  timer->startTimer("MPINSSolver2D - Init");
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
    zero_dat(dPdN_oi[0]);
    zero_dat(dPdN_oi[1]);
    zero_dat(n[0][0]);
    zero_dat(n[0][1]);
    zero_dat(n[1][0]);
    zero_dat(n[1][1]);
    zero_dat(pr);

    if(surface_tension) {
      zero_dat(st[0][0]);
      zero_dat(st[0][1]);
      zero_dat(st[1][0]);
      zero_dat(st[1][1]);
    }
  }

  sub_cycle_dt = h / (DG_ORDER * DG_ORDER * max_vel());
  if(!dt_forced) {
    dt = sub_cycle_dt;
    if(resuming && it_pre_sub_cycle <= 0 && sub_cycles > 1)
      dt = sub_cycle_dt * sub_cycles;
  } else {
    sub_cycle_dt = dt;
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

  lsSolver->init();
  lsSolver->set_bc_types(bc_types);
  lsSolver->getRhoMu(rho, mu);

  timer->endTimer("MPINSSolver2D - Init");
}

void MPINSSolver2D::step() {
  timer->startTimer("MPINSSolver2D - Advection");
  advection();
  timer->endTimer("MPINSSolver2D - Advection");

  timer->startTimer("MPINSSolver2D - Pressure");
  pressure();
  timer->endTimer("MPINSSolver2D - Pressure");

  timer->startTimer("MPINSSolver2D - Viscosity");
  viscosity();
  timer->endTimer("MPINSSolver2D - Viscosity");

  timer->startTimer("MPINSSolver2D - Surface");
  surface();
  timer->endTimer("MPINSSolver2D - Surface");

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

void MPINSSolver2D::surface_tension_grad(op_dat dx, op_dat dy) {
  DGTempDat heaviside = dg_dat_pool->requestTempDatCells(DG_NP);
  op_par_loop(ins_2d_st_new_0, "ins_2d_st_new_0", mesh->cells,
              op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(lsSolver->s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(heaviside.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->grad_over_int_with_central_flux(heaviside.dat, dx, dy);
  dg_dat_pool->releaseTempDatCells(heaviside);
}

void MPINSSolver2D::surface_tension_grad_over_int(op_dat dx, op_dat dy) {
  // grad heaviside
  DGTempDat st_tmp_0 = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, lsSolver->s, 0.0, st_tmp_0.dat);

  op_par_loop(ins_2d_st_0, "ins_2d_st_0", mesh->cells,
              op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(st_tmp_0.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, st_tmp_0.dat, 0.0, dx);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, st_tmp_0.dat, 0.0, dy);

  dg_dat_pool->releaseTempDatCells(st_tmp_0);

  op_par_loop(ins_2d_st_1, "ins_2d_st_1", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(dx, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dy, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  DGTempDat sM = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat sP = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(ins_2d_st_2, "ins_2d_st_2", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(lsSolver->s, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(sM.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(sP.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  if(mesh->bface2cells) {
    op_par_loop(ins_2d_st_3, "ins_2d_st_3", mesh->bfaces,
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(lsSolver->s, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(sM.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
                op_arg_dat(sP.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW));
  }

  DGTempDat sM_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  DGTempDat sP_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);

  timer->startTimer("INSSolverBase2D - advec_current_non_linear_over_int - Interp Surf");
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, sM.dat, 0.0, sM_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, sP.dat, 0.0, sP_cub.dat);

  dg_dat_pool->releaseTempDatCells(sM);
  dg_dat_pool->releaseTempDatCells(sP);
  timer->endTimer("INSSolverBase2D - advec_current_non_linear_over_int - Interp Surf");

  op_par_loop(ins_2d_st_4, "ins_2d_st_4", mesh->cells,
              op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nx_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(sM_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(sP_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, sM_cub.dat, -1.0, dx);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, sP_cub.dat, -1.0, dy);

  dg_dat_pool->releaseTempDatCells(sM_cub);
  dg_dat_pool->releaseTempDatCells(sP_cub);
}

void MPINSSolver2D::surface_tension_curvature(op_dat curv) {
  DGTempDat tmp_normal_x  = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_normal_y  = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->grad_over_int_with_central_flux(lsSolver->s, tmp_normal_x.dat, tmp_normal_y.dat);

  // Unit normals
  op_par_loop(ins_2d_st_5, "ins_2d_st_5", mesh->cells,
              op_arg_dat(tmp_normal_x.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_normal_y.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  mesh->div_over_int_with_central_flux(tmp_normal_x.dat, tmp_normal_y.dat, curv);

  // Curvature correction
  op_par_loop(ins_2d_st_8, "ins_2d_st_8", mesh->cells,
              op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(lsSolver->s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(curv, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  dg_dat_pool->releaseTempDatCells(tmp_normal_x);
  dg_dat_pool->releaseTempDatCells(tmp_normal_y);
}

// Calculate Nonlinear Terms
void MPINSSolver2D::advection() {
  if(surface_tension) {
    // Calculate surface tension
    // grad heaviside
    if(over_int_surface_tension)
      surface_tension_grad_over_int(st[currentInd][0], st[currentInd][1]);
    else
      surface_tension_grad(st[currentInd][0], st[currentInd][1]);

    // Calculate curvature
    DGTempDat tmp_curvature = dg_dat_pool->requestTempDatCells(DG_NP);
    surface_tension_curvature(tmp_curvature.dat);
/*
    // Apply curvature and weber number (placeholder currently)
    op_par_loop(ins_2d_st_6, "ins_2d_st_6", mesh->cells,
                op_arg_dat(tmp_curvature.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(st[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(st[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    dg_dat_pool->releaseTempDatCells(tmp_curvature);
*/
    DGTempDat rho_oi  = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
    lsSolver->getRhoVolOI(rho_oi.dat);
    DGTempDat curv_oi  = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, tmp_curvature.dat, 0.0, curv_oi.dat);
    dg_dat_pool->releaseTempDatCells(tmp_curvature);
    DGTempDat stx_oi = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
    DGTempDat sty_oi = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, st[currentInd][0], 0.0, stx_oi.dat);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, st[currentInd][1], 0.0, sty_oi.dat);
    op_par_loop(ins_2d_st_7, "ins_2d_st_7", mesh->cells,
                op_arg_dat(curv_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rho_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
                op_arg_dat(stx_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
                op_arg_dat(sty_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));

    op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PROJ, stx_oi.dat, 0.0, st[currentInd][0]);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PROJ, sty_oi.dat, 0.0, st[currentInd][1]);

    dg_dat_pool->releaseTempDatCells(curv_oi);
    dg_dat_pool->releaseTempDatCells(rho_oi);
    dg_dat_pool->releaseTempDatCells(stx_oi);
    dg_dat_pool->releaseTempDatCells(sty_oi);

  }

  if(time == 0.0 || sub_cycles < 1 || it_pre_sub_cycle != 0) {
    if(surface_tension)
      advec_standard(st[currentInd][0], st[currentInd][1], st[(currentInd + 1) % 2][0], st[(currentInd + 1) % 2][1]);
    else
      advec_standard();
  } else {
    if(surface_tension)
      advec_sub_cycle(st[currentInd][0], st[currentInd][1], st[(currentInd + 1) % 2][0], st[(currentInd + 1) % 2][1]);
    else
      advec_sub_cycle();
  }
}

void MPINSSolver2D::apply_pressure_neumann_bc(op_dat divVelT) {
  DGTempDat curlVel = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat gradCurlVel[2];
  gradCurlVel[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  gradCurlVel[1] = dg_dat_pool->requestTempDatCells(DG_NP);

  mesh->curl(vel[currentInd][0], vel[currentInd][1], curlVel.dat);

  op_par_loop(mp_ins_2d_pr_mu, "mp_ins_2d_pr_mu", mesh->cells,
              op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(curlVel.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  mesh->grad(curlVel.dat, gradCurlVel[0].dat, gradCurlVel[1].dat);

  dg_dat_pool->releaseTempDatCells(curlVel);

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
                op_arg_dat(dPdN[currentInd], 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));

  op_par_loop(mp_ins_pr_bc_dpdn_factors, "mp_ins_pr_bc_dpdn_factors", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT);
}

void MPINSSolver2D::apply_pressure_neumann_bc_oi(op_dat divVelT) {
  DGTempDat curlVel = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat gradCurlVel[2];
  gradCurlVel[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  gradCurlVel[1] = dg_dat_pool->requestTempDatCells(DG_NP);

  mesh->curl(vel[currentInd][0], vel[currentInd][1], curlVel.dat);

  DGTempDat mu_oi = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  lsSolver->getMuVolOI(mu_oi.dat);
  DGTempDat curlVel_oi = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, curlVel.dat, 0.0, curlVel_oi.dat);

  op_par_loop(mp_ins_2d_pr_mu_oi, "mp_ins_2d_pr_mu_oi", mesh->cells,
              op_arg_dat(mu_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(curlVel_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PROJ, curlVel_oi.dat, 0.0, curlVel.dat);
  dg_dat_pool->releaseTempDatCells(mu_oi);
  dg_dat_pool->releaseTempDatCells(curlVel_oi);

  mesh->grad(curlVel.dat, gradCurlVel[0].dat, gradCurlVel[1].dat);

  dg_dat_pool->releaseTempDatCells(curlVel);

  DGTempDat rho_oi = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  lsSolver->getRhoSurfOI(rho_oi.dat);
  op_par_loop(mp_ins_2d_pr_bc_oi_0, "mp_ins_2d_pr_bc_oi_0", mesh->bfaces,
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
              op_arg_dat(rho_oi.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN_oi[currentInd], 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC));

  dg_dat_pool->releaseTempDatCells(rho_oi);
  dg_dat_pool->releaseTempDatCells(gradCurlVel[0]);
  dg_dat_pool->releaseTempDatCells(gradCurlVel[1]);

  op_par_loop(mp_ins_2d_pr_bc_oi_1, "mp_ins_2d_pr_bc_oi_1", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN_oi[currentInd], -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN_oi[(currentInd + 1) % 2], -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, dPdN_oi[(currentInd + 1) % 2], 1.0, divVelT);
}

void MPINSSolver2D::update_pressure_matrices(DGTempDat &pr_factor) {
  pr_factor = dg_dat_pool->requestTempDatCells(DG_NP);
  op_par_loop(reciprocal, "reciprocal", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pr_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  FactorPoissonMatrixFreeDiag2D *tmpPressureMatrix = dynamic_cast<FactorPoissonMatrixFreeDiag2D*>(pressureMatrix);
  tmpPressureMatrix->set_factor(pr_factor.dat);
  pressureCoarseMatrix->set_factor(pr_factor.dat);
}

void MPINSSolver2D::update_pressure_matrices_oi(DGTempDat &pr_factor, DGTempDat &pr_factor_oi, DGTempDat &pr_factor_surf_oi) {
  pr_factor = dg_dat_pool->requestTempDatCells(DG_NP);
  pr_factor_oi = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  pr_factor_surf_oi = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);

  lsSolver->getRhoVolOI(pr_factor_oi.dat);
  op_par_loop(mp_ins_2d_pr_oi_0, "mp_ins_2d_pr_oi_0", mesh->cells,
              op_arg_dat(pr_factor_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));

  lsSolver->getRhoSurfOI(pr_factor_surf_oi.dat);
  op_par_loop(mp_ins_2d_pr_oi_1, "mp_ins_2d_pr_oi_1", mesh->cells,
              op_arg_dat(pr_factor_surf_oi.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PROJ, pr_factor_oi.dat, 0.0, pr_factor.dat);

  FactorPoissonMatrixFreeDiagOI2D *tmpPressureMatrix = dynamic_cast<FactorPoissonMatrixFreeDiagOI2D*>(pressureMatrix);
  tmpPressureMatrix->mat_free_set_factor_oi(pr_factor_oi.dat);
  tmpPressureMatrix->mat_free_set_factor_surf_oi(pr_factor_surf_oi.dat);
  tmpPressureMatrix->set_factor(pr_factor.dat);
  pressureCoarseMatrix->set_factor(pr_factor.dat);
}

void MPINSSolver2D::update_pressure_gradient(op_dat dpdx, op_dat dpdy) {
  op_par_loop(mp_ins_2d_pr_2, "mp_ins_2d_pr_2", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
}

void MPINSSolver2D::update_pressure_gradient_oi(op_dat dpdx, op_dat dpdy) {
  DGTempDat dpdx_oi = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat dpdy_oi = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, dpdx, 0.0, dpdx_oi.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, dpdy, 0.0, dpdy_oi.dat);

  DGTempDat rho_oi = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  lsSolver->getRhoVolOI(rho_oi.dat);

  op_par_loop(mp_ins_2d_pr_2_oi, "mp_ins_2d_pr_2_oi", mesh->cells,
              op_arg_dat(rho_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdx_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdy_oi.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));

  dg_dat_pool->releaseTempDatCells(rho_oi);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PROJ, dpdx_oi.dat, 0.0, dpdx);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PROJ, dpdy_oi.dat, 0.0, dpdy);

  dg_dat_pool->releaseTempDatCells(dpdx_oi);
  dg_dat_pool->releaseTempDatCells(dpdy_oi);
}

bool MPINSSolver2D::pressure() {
  timer->startTimer("MPINSSolver2D - Pressure RHS");
  DGTempDat divVelT = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->div_over_int_with_central_flux(velT[0], velT[1], divVelT.dat);

  // Calculate RHS of pressure solve
  const DG_FP div_factor = -1.0 / dt;
  op_par_loop(mp_ins_pr_div_factor, "mp_ins_pr_div_factor", mesh->cells,
              op_arg_gbl(&div_factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(divVelT.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  // Apply Neumann pressure boundary conditions
  if(mesh->bface2cells) {
    if(pr_over_int)
      apply_pressure_neumann_bc_oi(divVelT.dat);
    else
      apply_pressure_neumann_bc(divVelT.dat);
  }

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

  // Multiply RHS by mass matrix
  mesh->mass(divVelT.dat);
  timer->endTimer("MPINSSolver2D - Pressure RHS");

  // Call PETSc linear solver
  timer->startTimer("MPINSSolver2D - Pressure Linear Solve");
  DGTempDat pr_factor, pr_factor_oi, pr_factor_surf_oi;
  if(pr_over_int)
    update_pressure_matrices_oi(pr_factor_oi, pr_factor_surf_oi, pr_factor);
  else
    update_pressure_matrices(pr_factor);

  pressureMatrix->set_bc_types(pr_bc_types);
  pressureCoarseMatrix->set_bc_types(pr_bc_types);
  pressureSolver->set_coarse_matrix(pressureCoarseMatrix);
  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_bcs(bc_data);
  bool converged = pressureSolver->solve(divVelT.dat, pr);
  if(!converged)
    dg_abort("Pressure solve did not converge");

  if(pr_over_int) {
    dg_dat_pool->releaseTempDatCells(pr_factor_surf_oi);
    dg_dat_pool->releaseTempDatCells(pr_factor_oi);
  }
  dg_dat_pool->releaseTempDatCells(pr_factor);
  dg_dat_pool->releaseTempDatCells(divVelT);
  timer->endTimer("MPINSSolver2D - Pressure Linear Solve");

  if(pr_over_int)
    zero_dat(dPdN_oi[(currentInd + 1) % 2]);
  else
    zero_dat(dPdN[(currentInd + 1) % 2]);

  timer->startTimer("MPINSSolver2D - Pressure Projection");
  DGTempDat dpdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdy = dg_dat_pool->requestTempDatCells(DG_NP);

  // Calculate gradient of pressure
  mesh->grad_over_int_with_central_flux(pr, dpdx.dat, dpdy.dat);

  if(gravity_modified_pressure) {
    // Calculate average density
    DGTempDat tmp_rho = dg_dat_pool->requestTempDatCells(DG_NP);
    op_par_loop(copy_dg_np, "copy_dg_np", mesh->cells,
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_rho.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    mesh->mass(tmp_rho.dat);
    DG_FP avg_rho = 0.0;
    op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
                op_arg_gbl(&avg_rho, 1, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_rho.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
    avg_rho /= domain_area;
    dg_dat_pool->releaseTempDatCells(tmp_rho);

    // Add new gravity formulation
    op_par_loop(mp_ins_2d_pr_grav, "mp_ins_2d_pr_grav", mesh->cells,
                op_arg_gbl(&avg_rho, 1, DG_FP_STR, OP_READ),
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdx.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(dpdy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  }

  if(pr_over_int) {
    update_pressure_gradient_oi(dpdx.dat, dpdy.dat);
  } else {
    update_pressure_gradient(dpdx.dat, dpdy.dat);
  }

  project_velocity(dpdx.dat, dpdy.dat);

  dg_dat_pool->releaseTempDatCells(dpdx);
  dg_dat_pool->releaseTempDatCells(dpdy);
  timer->endTimer("MPINSSolver2D - Pressure Projection");

  return converged;
}

bool MPINSSolver2D::viscosity() {
  timer->startTimer("MPINSSolver2D - Viscosity RHS");
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

  DGTempDat vis_mm_factor = dg_dat_pool->requestTempDatCells(DG_NP);
  DG_FP factor  = reynolds / dt;
  DG_FP factor2 = g0 * reynolds / dt;
  op_par_loop(mp_ins_2d_vis_0, "mp_ins_2d_vis_0", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&factor2, 1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(visRHS[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(visRHS[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(vis_mm_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->mass(visRHS[0].dat);
  mesh->mass(visRHS[1].dat);
  timer->endTimer("MPINSSolver2D - Viscosity RHS");

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

    DGTempDat tmp_art_vis_x;
    if(shock_capturing) {
      tmp_art_vis_x = dg_dat_pool->requestTempDatCells(DG_NP);
      DGTempDat tmp_art_vis_y = dg_dat_pool->requestTempDatCells(DG_NP);
      calc_art_vis(velTT[0], tmp_art_vis_x.dat);
      calc_art_vis(velTT[1], tmp_art_vis_y.dat);

      op_par_loop(art_vis_2d_max, "art_vis_2d_max", mesh->cells,
                  op_arg_dat(tmp_art_vis_y.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_art_vis_x.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

      dg_dat_pool->releaseTempDatCells(tmp_art_vis_y);

      op_par_loop(mp_ins_2d_add_mu, "mp_ins_2d_add_mu", mesh->cells,
                  op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_art_vis_x.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
      
      slipViscousMatrix->set_factor(tmp_art_vis_x.dat);
    } else {
      slipViscousMatrix->set_factor(mu);
    }
    slipViscousMatrix->set_mm_factor(vis_mm_factor.dat);
    slipViscousMatrix->set_bc_types(vis_bc_types, vis_bc_types_2);
    slipViscousSolver->set_bcs(bc_data, bc_data_2);
    slipViscousSolver->set_inv_mass_recp_factor(vis_mm_factor.dat);
    slipViscousMatrix->calc_diag();
    slipViscousMatrix->calc_inv_block_diag();

    bool converged = slipViscousSolver->solve(visRHS[0].dat, visRHS[1].dat, vel[(currentInd + 1) % 2][0], vel[(currentInd + 1) % 2][1]);

    if(!converged)
      dg_abort("Viscosity solve did not converge");

    if(shock_capturing)
      dg_dat_pool->releaseTempDatCells(tmp_art_vis_x);
    dg_dat_pool->releaseTempDatCells(vis_mm_factor);
    dg_dat_pool->releaseTempDatCells(visRHS[0]);
    dg_dat_pool->releaseTempDatCells(visRHS[1]);

    return converged;
  } else {
    // Call PETSc linear solver
    timer->startTimer("MPINSSolver2D - Viscosity Linear Solve");
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

    DGTempDat tmp_art_vis = dg_dat_pool->requestTempDatCells(DG_NP);
    if(shock_capturing) {
      calc_art_vis(velTT[0], tmp_art_vis.dat);

      op_par_loop(mp_ins_2d_add_mu, "mp_ins_2d_add_mu", mesh->cells,
                  op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

      if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
        FactorMMPoissonMatrixFreeDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(tmp_art_vis.dat);
      } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
        FactorMMPoissonMatrixFreeBlockDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(tmp_art_vis.dat);
      }
    } else {
      if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
        FactorMMPoissonMatrixFreeDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(mu);
      } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
        FactorMMPoissonMatrixFreeBlockDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(mu);
      }
    }
    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag2D*>(viscosityMatrix);
      tmpMatrix->set_mm_factor(vis_mm_factor.dat);
      tmpMatrix->calc_mat_partial();
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag2D*>(viscosityMatrix);
      tmpMatrix->set_mm_factor(vis_mm_factor.dat);
      tmpMatrix->calc_mat_partial();
    }

    bool convergedX = viscositySolver->solve(visRHS[0].dat, vel[(currentInd + 1) % 2][0]);
    if(!convergedX)
      dg_abort("Viscosity X solve did not converge");

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

      op_par_loop(mp_ins_2d_add_mu, "mp_ins_2d_add_mu", mesh->cells,
                  op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

      if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
        FactorMMPoissonMatrixFreeDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(tmp_art_vis.dat);
      } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
        FactorMMPoissonMatrixFreeBlockDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(tmp_art_vis.dat);
      }
    } else {
      if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
        FactorMMPoissonMatrixFreeDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(mu);
      } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
        FactorMMPoissonMatrixFreeBlockDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag2D*>(viscosityMatrix);
        tmpMatrix->set_factor(mu);
      }
    }
    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag2D*>(viscosityMatrix);
      tmpMatrix->set_mm_factor(vis_mm_factor.dat);
      tmpMatrix->calc_mat_partial();
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag2D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag2D*>(viscosityMatrix);
      tmpMatrix->set_mm_factor(vis_mm_factor.dat);
      tmpMatrix->calc_mat_partial();
    }
    bool convergedY = viscositySolver->solve(visRHS[1].dat, vel[(currentInd + 1) % 2][1]);
    if(!convergedY)
      dg_abort("Viscosity Y solve did not converge");
    timer->endTimer("MPINSSolver2D - Viscosity Linear Solve");

    dg_dat_pool->releaseTempDatCells(tmp_art_vis);
    dg_dat_pool->releaseTempDatCells(vis_mm_factor);
    dg_dat_pool->releaseTempDatCells(visRHS[0]);
    dg_dat_pool->releaseTempDatCells(visRHS[1]);

    return convergedX && convergedY;
  }
}

void MPINSSolver2D::dump_checkpoint_data(const std::string &filename) {
  INSSolverBase2D::dump_checkpoint_data(filename);

  if(surface_tension) {
    op_fetch_data_hdf5_file(st[0][0], filename.c_str());
    op_fetch_data_hdf5_file(st[0][1], filename.c_str());
    op_fetch_data_hdf5_file(st[1][0], filename.c_str());
    op_fetch_data_hdf5_file(st[1][1], filename.c_str());
  }
  op_fetch_data_hdf5_file(lsSolver->s, filename.c_str());

  // TODO save constants in same HDF5 file
}

void MPINSSolver2D::dump_visualisation_data(const std::string &filename) {
  INSSolverBase2D::dump_visualisation_data(filename);

  if(values_to_save.count("surface_tension") != 0 && surface_tension) {
    op_fetch_data_hdf5_file(st[(currentInd + 1) % 2][0], filename.c_str());
    op_fetch_data_hdf5_file(st[(currentInd + 1) % 2][1], filename.c_str());
  }

  if(values_to_save.count("mu") != 0) {
    if(shock_capturing) {
      DGTempDat tmp_art_vis = dg_dat_pool->requestTempDatCells(DG_NP);
      calc_art_vis(velTT[0], tmp_art_vis.dat);

      op_par_loop(mp_ins_2d_add_mu, "copy_dg_np", mesh->cells,
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
    op_fetch_data_hdf5_file(lsSolver->kink, filename.c_str());
  }

  if(values_to_save.count("curvature") != 0) {
    op_dat curv = op_decl_dat_temp(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_curvature");
    surface_tension_curvature(curv);
    op_fetch_data_hdf5_file(curv, filename.c_str());
    op_free_dat_temp(curv);
  }
}

void MPINSSolver2D::surface() {
  lsSolver->setVelField(vel[(currentInd + 1) % 2][0], vel[(currentInd + 1) % 2][1]);
  if(it_pre_sub_cycle > 1 || sub_cycles < 1)
    lsSolver->step(dt, 1);
  else
    lsSolver->step(sub_cycle_dt, sub_cycles);
  lsSolver->getRhoMu(rho, mu);
}

op_dat MPINSSolver2D::get_ls() {
  return lsSolver->s;
}

DG_FP MPINSSolver2D::get_ls_alpha() {
  return lsSolver->alpha;
}
