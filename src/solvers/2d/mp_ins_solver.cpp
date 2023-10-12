#include "solvers/2d/mp_ins_solver.h"

// Include OP2 stuff
#include "op_seq.h"

#include <iostream>
#include <limits>

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"

#include "timing.h"
#include "config.h"

extern Timing *timer;
extern Config *config;
extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;

using namespace std;

MPINSSolver2D::MPINSSolver2D(DGMesh2D *m) : INSSolverBase2D(m) {
  resuming = false;

  setup_common();

  lsSolver = new LevelSetSolver2D(m);

  currentInd = 0;

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;
}

MPINSSolver2D::MPINSSolver2D(DGMesh2D *m, const std::string &filename, const int iter) : INSSolverBase2D(m, filename) {
  resuming = true;

  setup_common();

  lsSolver = new LevelSetSolver2D(mesh, filename);

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

void MPINSSolver2D::setup_common() {
  int tmp_st = 0;
  config->getInt("solver-options", "surface_tension", tmp_st);
  surface_tension = tmp_st == 1;
  double tmp_dt = -1.0;
  config->getDouble("solver-options", "force_dt", tmp_dt);
  dt_forced = tmp_dt > 0.0;
  if(dt_forced) dt = tmp_dt;

  pressureMatrix = new FactorPoissonMatrixFreeDiag2D(mesh);
  pressureCoarseMatrix = new FactorPoissonCoarseMatrix2D(mesh);
  viscosityMatrix = new FactorMMPoissonMatrixFreeDiag2D(mesh);
  pressureSolver = new PETScPMultigrid(mesh);
  pressureSolver->set_coarse_matrix(pressureCoarseMatrix);
  viscositySolver = new PETScJacobiSolver(mesh);

  int pr_tmp = 0;
  int vis_tmp = 0;
  config->getInt("solver-options", "pr_nullspace", pr_tmp);
  config->getInt("solver-options", "vis_nullspace", vis_tmp);

  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_nullspace(pr_tmp == 1);
  viscositySolver->set_matrix(viscosityMatrix);
  viscositySolver->set_nullspace(vis_tmp == 1);

  rho = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_rho");
  mu  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_mu");

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
  delete viscosityMatrix;
  delete pressureSolver;
  delete viscositySolver;
}

void MPINSSolver2D::init(const DG_FP re, const DG_FP refVel) {
  timer->startTimer("MPINSSolver2D - Init");
  INSSolverBase2D::init(re, refVel);

  reynolds = re;

  lsSolver->init();

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

    op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                op_arg_dat(n[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                op_arg_dat(n[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                op_arg_dat(n[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                op_arg_dat(n[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

    if(surface_tension) {
      op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                  op_arg_dat(st[0][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
      op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                  op_arg_dat(st[0][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
      op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                  op_arg_dat(st[1][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
      op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                  op_arg_dat(st[1][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    }

    op_par_loop(zero_np_1, "zero_np_1", mesh->cells,
                op_arg_dat(pr, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
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
                op_arg_dat(mesh->node_coords, -3, mesh->bface2nodes, 3, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,     -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(pr_bc_types,  -1, OP_ID, 1, "int", OP_WRITE),
                op_arg_dat(vis_bc_types, -1, OP_ID, 1, "int", OP_WRITE));
  }

  pressureCoarseMatrix->set_bc_types(pr_bc_types);
  pressureMatrix->set_bc_types(pr_bc_types);
  pressureSolver->init();
  viscositySolver->init();

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
  time += dt;

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
void MPINSSolver2D::advection() {
  if(surface_tension) {
    // Calculate surface tension
    // grad heaviside
    DGTempDat st_tmp_0 = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, lsSolver->s, 0.0, st_tmp_0.dat);

    op_par_loop(ins_2d_st_0, "ins_2d_st_0", mesh->cells,
                op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
                op_arg_dat(st_tmp_0.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));
    
    op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, st_tmp_0.dat, 0.0, st[currentInd][0]);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, st_tmp_0.dat, 0.0, st[currentInd][1]);

    dg_dat_pool->releaseTempDatCells(st_tmp_0);

    op_par_loop(ins_2d_st_1, "ins_2d_st_1", mesh->cells,
                op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
                op_arg_dat(st[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(st[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

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

    op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, sM_cub.dat, -1.0, st[currentInd][0]);
    op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_LIFT, sP_cub.dat, -1.0, st[currentInd][1]);

    dg_dat_pool->releaseTempDatCells(sM_cub);
    dg_dat_pool->releaseTempDatCells(sP_cub);

    // Calculate curvature
    DGTempDat tmp_normal_x  = dg_dat_pool->requestTempDatCells(DG_NP);
    DGTempDat tmp_normal_y  = dg_dat_pool->requestTempDatCells(DG_NP);
    DGTempDat tmp_curvature = dg_dat_pool->requestTempDatCells(DG_NP);
    mesh->grad_with_central_flux(lsSolver->s, tmp_normal_x.dat, tmp_normal_y.dat);

    // Unit normals
    op_par_loop(ins_2d_st_5, "ins_2d_st_5", mesh->cells,
                op_arg_dat(tmp_normal_x.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(tmp_normal_y.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    mesh->div_with_central_flux(tmp_normal_x.dat, tmp_normal_y.dat, tmp_curvature.dat);

    // Apply curvature and weber number (placeholder currently)
    op_par_loop(ins_2d_st_6, "ins_2d_st_6", mesh->cells,
                op_arg_gbl(constants->decrease_order_ptr, DG_NP * DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_curvature.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(st[currentInd][0], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
                op_arg_dat(st[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
    
    dg_dat_pool->releaseTempDatCells(tmp_normal_x);
    dg_dat_pool->releaseTempDatCells(tmp_normal_y);
    dg_dat_pool->releaseTempDatCells(tmp_curvature);
  }

  if(time == 0.0 || sub_cycles < 1 || it_pre_sub_cycle != 0) {
    if(surface_tension)
      advec_standard(st[currentInd][0], st[currentInd][1], st[(currentInd + 1) % 2][0], st[(currentInd + 1) % 2][1]);
    else
      advec_standard();
  } else {
    advec_sub_cycle();
  }
}

bool MPINSSolver2D::pressure() {
  timer->startTimer("MPINSSolver2D - Pressure RHS");
  DGTempDat divVelT = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat curlVel = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat gradCurlVel[2];
  gradCurlVel[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  gradCurlVel[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->div_with_central_flux(velT[0], velT[1], divVelT.dat);
  mesh->curl(vel[currentInd][0], vel[currentInd][1], curlVel.dat);

  op_par_loop(mp_ins_2d_pr_mu, "mp_ins_3d_pr_mu", mesh->cells,
              op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(curlVel.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

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

  DGTempDat prBC = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  // Apply Dirichlet BCs
  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(prBC.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  // Calculate RHS of pressure solve
  DGTempDat pr_factor = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(mp_ins_2d_pr_1, "mp_ins_2d_pr_1", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&dt, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
              op_arg_dat(divVelT.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pr_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT.dat);
  mesh->mass(divVelT.dat);
  timer->endTimer("MPINSSolver2D - Pressure RHS");

  // Call PETSc linear solver
  timer->startTimer("MPINSSolver2D - Pressure Linear Solve");
  pressureMatrix->set_factor(pr_factor.dat);
  pressureCoarseMatrix->set_factor(pr_factor.dat);
  pressureMatrix->set_bc_types(pr_bc_types);
  pressureCoarseMatrix->set_bc_types(pr_bc_types);
  pressureSolver->set_coarse_matrix(pressureCoarseMatrix);
  pressureSolver->set_matrix(pressureMatrix);
  pressureSolver->set_bcs(prBC.dat);
  bool converged = pressureSolver->solve(divVelT.dat, pr);

  dg_dat_pool->releaseTempDatCells(pr_factor);
  dg_dat_pool->releaseTempDatCells(divVelT);
  dg_dat_pool->releaseTempDatCells(prBC);
  timer->endTimer("MPINSSolver2D - Pressure Linear Solve");

  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  timer->startTimer("MPINSSolver2D - Pressure Projection");
  DGTempDat dpdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdy = dg_dat_pool->requestTempDatCells(DG_NP);

  // Calculate gradient of pressure
  mesh->grad_over_int_with_central_flux(pr, dpdx.dat, dpdy.dat);

  op_par_loop(mp_ins_2d_pr_2, "mp_ins_2d_pr_2", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdx.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  project_velocity(dpdx.dat, dpdy.dat);

  if(surface_tension) {
    filter(velTT[0]);
    filter(velTT[1]);
  }

  dg_dat_pool->releaseTempDatCells(dpdx);
  dg_dat_pool->releaseTempDatCells(dpdy);
  timer->endTimer("MPINSSolver2D - Pressure Projection");

  return converged;
}

bool MPINSSolver2D::viscosity() {
  timer->startTimer("MPINSSolver2D - Viscosity RHS");
  DG_FP time_n1 = time + dt;

  DGTempDat visBC[2];
  visBC[0] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  visBC[1] = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_2, "zero_npf_2", mesh->cells,
              op_arg_dat(visBC[0].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(visBC[1].dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  // Get BCs for viscosity solve
  if(mesh->bface2cells) {
    op_par_loop(ins_vis_bc_2d, "ins_vis_bc_2d", mesh->bfaces,
                op_arg_gbl(&time_n1, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->x,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(visBC[0].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC),
                op_arg_dat(visBC[1].dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

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

  // Call PETSc linear solver
  timer->startTimer("MPINSSolver2D - Viscosity Linear Solve");
  viscosityMatrix->set_factor(mu);
  viscosityMatrix->set_mm_factor(vis_mm_factor.dat);
  viscosityMatrix->set_bc_types(vis_bc_types);
  viscosityMatrix->calc_mat_partial();
  viscositySolver->set_bcs(visBC[0].dat);
  bool convergedX = viscositySolver->solve(visRHS[0].dat, vel[(currentInd + 1) % 2][0]);

  viscositySolver->set_bcs(visBC[1].dat);
  bool convergedY = viscositySolver->solve(visRHS[1].dat, vel[(currentInd + 1) % 2][1]);
  timer->endTimer("MPINSSolver2D - Viscosity Linear Solve");

  dg_dat_pool->releaseTempDatCells(vis_mm_factor);
  dg_dat_pool->releaseTempDatCells(visRHS[0]);
  dg_dat_pool->releaseTempDatCells(visRHS[1]);
  dg_dat_pool->releaseTempDatCells(visBC[0]);
  dg_dat_pool->releaseTempDatCells(visBC[1]);

  return convergedX && convergedY;
}

void MPINSSolver2D::dump_data(const std::string &filename) {
  timer->startTimer("MPINSSolver2D - Dump Data");
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
  
  if(surface_tension) {
    op_fetch_data_hdf5_file(st[0][0], filename.c_str());
    op_fetch_data_hdf5_file(st[0][1], filename.c_str());
    op_fetch_data_hdf5_file(st[1][0], filename.c_str());
    op_fetch_data_hdf5_file(st[1][1], filename.c_str());
  }
  op_fetch_data_hdf5_file(dPdN[0], filename.c_str());
  op_fetch_data_hdf5_file(dPdN[1], filename.c_str());
  op_fetch_data_hdf5_file(velT[0], filename.c_str());
  op_fetch_data_hdf5_file(velT[1], filename.c_str());
  // mesh->grad_over_int_with_central_flux(pr, velTT[0], velTT[1]);
  op_fetch_data_hdf5_file(velTT[0], filename.c_str());
  op_fetch_data_hdf5_file(velTT[1], filename.c_str());
  op_fetch_data_hdf5_file(pr, filename.c_str());
  op_fetch_data_hdf5_file(mu, filename.c_str());
  op_fetch_data_hdf5_file(rho, filename.c_str());
  op_fetch_data_hdf5_file(lsSolver->s, filename.c_str());
  timer->endTimer("MPINSSolver2D - Dump Data");
}

void MPINSSolver2D::surface() {
  lsSolver->setVelField(vel[(currentInd + 1) % 2][0], vel[(currentInd + 1) % 2][1]);
  lsSolver->step(dt);
  lsSolver->getRhoMu(rho, mu);
}
