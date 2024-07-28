#include "solvers/3d/mp_ins_solver.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_utils.h"
#include "dg_constants/dg_constants.h"

#include "timing.h"
#include "config.h"
#include "op2_utils.h"
#include "dg_matrices/3d/factor_mm_poisson_matrix_free_diag_3d.h"
#include "dg_matrices/3d/factor_mm_poisson_matrix_free_block_diag_3d.h"
#include "dg_linear_solvers/petsc_jacobi.h"
#include "dg_linear_solvers/petsc_block_jacobi.h"
#include "dg_linear_solvers/initial_guess_extrapolation.h"
#include "dg_abort.h"

#include <string>

extern Timing *timer;
extern DGConstants *constants;
extern Config *config;
extern DGDatPool *dg_dat_pool;

MPINSSolver3D::MPINSSolver3D(DGMesh3D *m, const DG_FP re) : INSSolverBase3D(m) {
  resuming = false;

  a0 = 1.0;
  a1 = 0.0;
  b0 = 1.0;
  b1 = 0.0;
  g0 = 1.0;

  dt = 0.0;

  reynolds = re;

  lsSolver = new LevelSetSolver3D(mesh);

  setup_common();
}

MPINSSolver3D::MPINSSolver3D(DGMesh3D *m, const DG_FP re, const std::string &filename, const int iter) : INSSolverBase3D(m, filename) {
  resuming = true;

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

  lsSolver = new LevelSetSolver3D(mesh, filename);

  setup_common();
}

void MPINSSolver3D::setup_common() {
  int tmp_st = 0;
  config->getInt("solver-options", "surface_tension", tmp_st);
  surface_tension = tmp_st == 1;
  int tmp_grav = 0;
  config->getInt("solver-options", "gravity_modified_pressure", tmp_grav);
  gravity_modified_pressure = tmp_grav == 1;
  int tmp_force_vel = 0;
  config->getInt("force-superficial-velocity", "on", tmp_force_vel);
  force_superficial_velocity = tmp_force_vel == 1;
  fsv_relaxation_factor = 0.9;
  config->getDouble("force-superficial-velocity", "relaxation_factor", fsv_relaxation_factor);
  fsv_factor = 750.0;
  config->getDouble("force-superficial-velocity", "initial_forcing_value", fsv_factor);

  if(gravity && gravity_modified_pressure)
    dg_abort("Do not use both \'gravity\' and \'gravity_modified_pressure\'");

  // Pressure matrix and solver
  std::string pr_solver = "p-multigrid";
  config->getStr("pressure-solve", "preconditioner", pr_solver);
  pressureSolverType = set_solver_type(pr_solver);
  if(pressureSolverType != LinearSolver::PETSC_PMULTIGRID)
    dg_abort("Only \'p-multigrid\' preconditioner is supported for 3D multiphase flow.");
  int tmp_pr_over_int = 0;
  config->getInt("pressure-solve", "over_int", tmp_pr_over_int);
  pr_over_int = tmp_pr_over_int != 0;
  if(pr_over_int) {
    pressureMatrix = new FactorPoissonMatrixFreeDiagOI3D(mesh);
  } else {
    pressureMatrix = new FactorPoissonMatrixFreeDiag3D(mesh);
  }
  coarsePressureMatrix = new FactorPoissonCoarseMatrix3D(mesh);
  pressureSolver = new PETScPMultigrid(mesh);

  // Viscous matrix and solver
  std::string vis_solver = "jacobi";
  config->getStr("viscous-solve", "preconditioner", vis_solver);
  viscositySolverType = set_solver_type(vis_solver);
  if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
    viscosityMatrix = new FactorMMPoissonMatrixFreeDiag3D(mesh);
    viscositySolver = new PETScJacobiSolver(mesh);
  } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
    viscosityMatrix = new FactorMMPoissonMatrixFreeBlockDiag3D(mesh);
    viscositySolver = new PETScBlockJacobiSolver(mesh);
  } else {
    dg_abort("Only \'jacobi\' and \'block-jacobi\' preconditioner is supported for 3D multiphase flow.");
  }

  pressureSolver->set_matrix(pressureMatrix);
  viscositySolver->set_matrix(viscosityMatrix);

  setup_pressure_viscous_solvers(pressureSolver, viscositySolver);

  rho = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_rho");
  mu  = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_mu");

  tmp_npf_bc = op_decl_dat(mesh->bfaces, DG_NPF, DG_FP_STR, (DG_FP *)NULL, "ins_solver_tmp_npf_bc");
  tmp_bc_1 = op_decl_dat(mesh->bfaces, 1, "int", (int *)NULL, "ins_solver_tmp_bc_1");
  art_vis  = op_decl_dat(mesh->cells, 1, DG_FP_STR, (DG_FP *)NULL, "ins_solver_art_vis");

  if(surface_tension) {
    st[0][0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st00");
    st[0][1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st01");
    st[0][2] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st02");
    st[1][0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st10");
    st[1][1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st11");
    st[1][2] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ins_solver_st12");
  }

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

void MPINSSolver3D::init() {
  timer->startTimer("MPINSSolver3D - Init");
  INSSolverBase3D::init();

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
  if(resuming && it_pre_sub_cycle <= 0 && sub_cycles > 1)
    dt = sub_cycle_dt * sub_cycles;
  // dt *= 1e-2;
  op_printf("INS dt is %g\n", dt);

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

    if(surface_tension) {
      zero_dat(st[0][0]);
      zero_dat(st[0][1]);
      zero_dat(st[0][2]);
      zero_dat(st[1][0]);
      zero_dat(st[1][1]);
      zero_dat(st[1][2]);
    }
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

void MPINSSolver3D::surface_tension_grad(op_dat dx, op_dat dy, op_dat dz) {
  // grad heaviside
  DGTempDat st_tmp_0 = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, lsSolver->s, 0.0, st_tmp_0.dat);

  op_par_loop(ins_3d_st_0, "ins_3d_st_0", mesh->cells,
              op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(st_tmp_0.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PDR, st_tmp_0.dat, 0.0, dx);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PDS, st_tmp_0.dat, 0.0, dy);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PDT, st_tmp_0.dat, 0.0, dz);

  dg_dat_pool->releaseTempDatCells(st_tmp_0);

  op_par_loop(ins_3d_st_1, "ins_3d_st_1", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 10, DG_FP_STR, OP_READ),
              op_arg_dat(dx, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dy, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dz, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  DGTempDat fX_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_3D_NP);
  DGTempDat fY_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_3D_NP);
  DGTempDat fZ_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_3D_NP);

  zero_dat(fX_cub.dat);
  zero_dat(fY_cub.dat);
  zero_dat(fZ_cub.dat);

  op_par_loop(ins_3d_st_2, "ins_3d_st_2", mesh->faces,
              op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(lsSolver->s, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(fX_cub.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(fY_cub.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(fZ_cub.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_WRITE));

  if(mesh->bface2cells) {
    op_par_loop(ins_3d_st_3, "ins_3d_st_3", mesh->bfaces,
                op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(lsSolver->s, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(fX_cub.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_INC),
                op_arg_dat(fY_cub.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_INC),
                op_arg_dat(fZ_cub.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF3D_LIFT, fX_cub.dat, -1.0, dx);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF3D_LIFT, fY_cub.dat, -1.0, dy);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF3D_LIFT, fZ_cub.dat, -1.0, dz);

  dg_dat_pool->releaseTempDatCells(fX_cub);
  dg_dat_pool->releaseTempDatCells(fY_cub);
  dg_dat_pool->releaseTempDatCells(fZ_cub);
}

void MPINSSolver3D::surface_tension_curvature(op_dat curv) {
  DGTempDat tmp_normal_x  = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_normal_y  = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_normal_z  = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->grad_with_central_flux(lsSolver->s, tmp_normal_x.dat, tmp_normal_y.dat, tmp_normal_z.dat);

  // Unit normals
  op_par_loop(ins_3d_st_4, "ins_3d_st_4", mesh->cells,
              op_arg_dat(tmp_normal_x.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_normal_y.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_normal_z.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  mesh->div_with_central_flux(tmp_normal_x.dat, tmp_normal_y.dat, tmp_normal_z.dat, curv);

  // TODO Curvature correction

  dg_dat_pool->releaseTempDatCells(tmp_normal_x);
  dg_dat_pool->releaseTempDatCells(tmp_normal_y);
  dg_dat_pool->releaseTempDatCells(tmp_normal_z);
}

void MPINSSolver3D::advection() {
  if(surface_tension) {
    // Calculate surface tension
    // grad heaviside
    surface_tension_grad(st[currentInd][0], st[currentInd][1], st[currentInd][2]);

    // Calculate curvature
    DGTempDat tmp_curvature = dg_dat_pool->requestTempDatCells(DG_NP);
    surface_tension_curvature(tmp_curvature.dat);
    DGTempDat rho_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
    lsSolver->getRhoVolOI(rho_oi.dat);
    DGTempDat curv_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, tmp_curvature.dat, 0.0, curv_oi.dat);
    dg_dat_pool->releaseTempDatCells(tmp_curvature);
    DGTempDat stx_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
    DGTempDat sty_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
    DGTempDat stz_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, st[currentInd][0], 0.0, stx_oi.dat);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, st[currentInd][1], 0.0, sty_oi.dat);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, st[currentInd][2], 0.0, stz_oi.dat);

    op_par_loop(ins_3d_st_5, "ins_3d_st_5", mesh->cells,
                op_arg_dat(curv_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rho_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
                op_arg_dat(stx_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW),
                op_arg_dat(sty_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW),
                op_arg_dat(stz_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW));

    op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PROJ, stx_oi.dat, 0.0, st[currentInd][0]);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PROJ, sty_oi.dat, 0.0, st[currentInd][1]);
    op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PROJ, stz_oi.dat, 0.0, st[currentInd][2]);

    dg_dat_pool->releaseTempDatCells(curv_oi);
    dg_dat_pool->releaseTempDatCells(rho_oi);
    dg_dat_pool->releaseTempDatCells(stx_oi);
    dg_dat_pool->releaseTempDatCells(sty_oi);
    dg_dat_pool->releaseTempDatCells(stz_oi);
  }

  if(time == 0.0 || sub_cycles < 1 || it_pre_sub_cycle != 0) {
    if(surface_tension)
      advec_standard(st[currentInd][0], st[currentInd][1], st[currentInd][2],
                     st[(currentInd + 1) % 2][0], st[(currentInd + 1) % 2][1],
                     st[(currentInd + 1) % 2][2]);
    else
      advec_standard();
  } else {
    if(surface_tension)
      advec_sub_cycle(st[currentInd][0], st[currentInd][1], st[currentInd][2],
                      st[(currentInd + 1) % 2][0], st[(currentInd + 1) % 2][1],
                      st[(currentInd + 1) % 2][2]);
    else
      advec_sub_cycle();
  }
}

void MPINSSolver3D::apply_pressure_neumann_bc(op_dat divVelT) {
  DGTempDat curlVel[3];
  curlVel[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  curlVel[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  curlVel[2] = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat curl2Vel[3];
  curl2Vel[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  curl2Vel[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  curl2Vel[2] = dg_dat_pool->requestTempDatCells(DG_NP);
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

  op_par_loop(mp_ins_pr_bc_dpdn_factors, "mp_ins_pr_bc_dpdn_factors", mesh->cells,
              op_arg_gbl(&b0, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&b1, 1, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[currentInd], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_READ),
              op_arg_dat(dPdN[(currentInd + 1) % 2], -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, dPdN[(currentInd + 1) % 2], 1.0, divVelT);

  dg_dat_pool->releaseTempDatCells(curl2Vel[0]);
  dg_dat_pool->releaseTempDatCells(curl2Vel[1]);
  dg_dat_pool->releaseTempDatCells(curl2Vel[2]);
}

void MPINSSolver3D::apply_pressure_neumann_bc_oi(op_dat divVelT) {
  // Temp
  apply_pressure_neumann_bc(divVelT);
}

void MPINSSolver3D::update_pressure_matrices(DGTempDat &pr_factor) {
  pr_factor = dg_dat_pool->requestTempDatCells(DG_NP);
  op_par_loop(reciprocal, "reciprocal", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(pr_factor.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  FactorPoissonMatrixFreeDiag3D *tmpPressureMatrix = dynamic_cast<FactorPoissonMatrixFreeDiag3D*>(pressureMatrix);
  tmpPressureMatrix->set_factor(pr_factor.dat);
  coarsePressureMatrix->set_factor(pr_factor.dat);
}

void MPINSSolver3D::update_pressure_matrices_oi(DGTempDat &pr_factor, DGTempDat &pr_factor_oi, DGTempDat &pr_factor_surf_oi) {
  pr_factor = dg_dat_pool->requestTempDatCells(DG_NP);
  pr_factor_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  pr_factor_surf_oi = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_3D_NP);

  lsSolver->getRhoVolOI(pr_factor_oi.dat);
  op_par_loop(mp_ins_3d_pr_fact_oi_0, "mp_ins_3d_pr_fact_oi_0", mesh->cells,
              op_arg_dat(pr_factor_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW));

  lsSolver->getRhoSurfOI(pr_factor_surf_oi.dat);
  op_par_loop(mp_ins_3d_pr_fact_oi_1, "mp_ins_3d_pr_fact_oi_1", mesh->cells,
              op_arg_dat(pr_factor_surf_oi.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PROJ, pr_factor_oi.dat, 0.0, pr_factor.dat);

  FactorPoissonMatrixFreeDiagOI3D *tmpPressureMatrix = dynamic_cast<FactorPoissonMatrixFreeDiagOI3D*>(pressureMatrix);
  tmpPressureMatrix->mat_free_set_factor_oi(pr_factor_oi.dat);
  tmpPressureMatrix->mat_free_set_factor_surf_oi(pr_factor_surf_oi.dat);
  tmpPressureMatrix->set_factor(pr_factor.dat);
  coarsePressureMatrix->set_factor(pr_factor.dat);
}

void MPINSSolver3D::update_pressure_gradient(op_dat dpdx, op_dat dpdy, op_dat dpdz) {
  op_par_loop(mp_ins_3d_pr_2, "mp_ins_3d_pr_2", mesh->cells,
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdx, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdy, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdz, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
}

void MPINSSolver3D::update_pressure_gradient_oi(op_dat dpdx, op_dat dpdy, op_dat dpdz) {
  DGTempDat dpdx_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  DGTempDat dpdy_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  DGTempDat dpdz_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, dpdx, 0.0, dpdx_oi.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, dpdy, 0.0, dpdy_oi.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, dpdz, 0.0, dpdz_oi.dat);

  DGTempDat rho_oi = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  lsSolver->getRhoVolOI(rho_oi.dat);

  op_par_loop(mp_ins_3d_pr_2_oi, "mp_ins_3d_pr_2_oi", mesh->cells,
              op_arg_dat(rho_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dpdx_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdy_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(dpdz_oi.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PROJ, dpdx_oi.dat, 0.0, dpdx);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PROJ, dpdy_oi.dat, 0.0, dpdy);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PROJ, dpdz_oi.dat, 0.0, dpdz);

  dg_dat_pool->releaseTempDatCells(rho_oi);
  dg_dat_pool->releaseTempDatCells(dpdx_oi);
  dg_dat_pool->releaseTempDatCells(dpdy_oi);
  dg_dat_pool->releaseTempDatCells(dpdz_oi);
}

void MPINSSolver3D::pressure() {
  DGTempDat divVelT = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->div_with_central_flux(velT[0], velT[1], velT[2], divVelT.dat);

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

  // Multiply RHS by mass matrix
  mesh->mass(divVelT.dat);

  // Call PETSc linear solver
  timer->startTimer("MPINSSolver3D - Pressure Linear Solve");
  DGTempDat pr_factor, pr_factor_oi, pr_factor_surf_oi;
  if(pr_over_int)
    update_pressure_matrices_oi(pr_factor_oi, pr_factor_surf_oi, pr_factor);
  else
    update_pressure_matrices(pr_factor);

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

  if(pr_over_int) {
    dg_dat_pool->releaseTempDatCells(pr_factor_surf_oi);
    dg_dat_pool->releaseTempDatCells(pr_factor_oi);
  }
  dg_dat_pool->releaseTempDatCells(pr_factor);
  dg_dat_pool->releaseTempDatCells(divVelT);
  timer->endTimer("MPINSSolver3D - Pressure Linear Solve");

  zero_dat(dPdN[(currentInd + 1) % 2]);

  timer->startTimer("MPINSSolver3D - Pressure Projection");
  DGTempDat dpdx = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdy = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat dpdz = dg_dat_pool->requestTempDatCells(DG_NP);
  mesh->grad_with_central_flux(pr, dpdx.dat, dpdy.dat, dpdz.dat);

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
    op_par_loop(mp_ins_pr_grav, "mp_ins_pr_grav", mesh->cells,
                op_arg_gbl(&avg_rho, 1, DG_FP_STR, OP_READ),
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  }

  // Currently only in the case where flow direction is (0,1)
  if(force_superficial_velocity) {
    // Volume that one phase occupies
    DGTempDat tmp_vol = dg_dat_pool->requestTempDatCells(DG_NP);
    op_par_loop(mp_ins_vol_phase, "mp_ins_vol_phase", mesh->cells,
                op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
                op_arg_dat(lsSolver->s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_vol.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    mesh->mass(tmp_vol.dat);
    DG_FP volume = 0.0;
    op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
                op_arg_gbl(&volume, 1, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_vol.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
    dg_dat_pool->releaseTempDatCells(tmp_vol);

    // Calculate the superficial velocity of one phase (sort of,
    // over the entire volume instead of a cross-section)
    DGTempDat tmp_sv = dg_dat_pool->requestTempDatCells(DG_NP);
    op_par_loop(mp_ins_sv_phase, "mp_ins_sv_phase", mesh->cells,
                op_arg_gbl(&lsSolver->alpha, 1, DG_FP_STR, OP_READ),
                op_arg_dat(lsSolver->s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(vel[currentInd][1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_sv.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
    mesh->mass(tmp_sv.dat);
    DG_FP sv = 0.0;
    op_par_loop(sum_dg_np, "sum_dg_np", mesh->cells,
                op_arg_gbl(&sv, 1, DG_FP_STR, OP_INC),
                op_arg_dat(tmp_sv.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
    dg_dat_pool->releaseTempDatCells(tmp_sv);
    sv /= volume;

    // Update forcing term
    fsv_factor -= fsv_relaxation_factor * (1.0 - sv);

    // Add forcing term
    // TODO account for g0 here??

    op_par_loop(mp_ins_pr_fsv, "mp_ins_pr_fsv", mesh->cells,
                op_arg_gbl(&fsv_factor, 1, DG_FP_STR, OP_READ),
                op_arg_dat(dpdy.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    op_printf("Air volume: %g\nAir superficial velocity: %g\nForcing factor: %g\n", volume, sv, fsv_factor);
  }

  if(pr_over_int) {
    update_pressure_gradient_oi(dpdx.dat, dpdy.dat, dpdz.dat);
  } else {
    update_pressure_gradient(dpdx.dat, dpdy.dat, dpdz.dat);
  }

  project_velocity(dpdx.dat, dpdy.dat, dpdz.dat);

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
  viscosityMatrix->set_bc_types(vis_bc_types);
  if(shock_capturing) {
    calc_art_vis(velTT[0], tmp_art_vis.dat);

    op_par_loop(mp_ins_3d_add_mu, "mp_ins_3d_add_mu", mesh->cells,
                op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
    }
  } else {
    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(mu);
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(mu);
    }
  }
  if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
    FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
    tmpMatrix->set_mm_factor(vis_mm_factor.dat);
    tmpMatrix->calc_mat_partial();
  } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
    FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
    tmpMatrix->set_mm_factor(vis_mm_factor.dat);
    tmpMatrix->calc_mat_partial();
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

  viscositySolver->set_bcs(vis_bc);
  if(shock_capturing) {
    calc_art_vis(velTT[1], tmp_art_vis.dat);

    op_par_loop(mp_ins_3d_add_mu, "mp_ins_3d_add_mu", mesh->cells,
                op_arg_dat(mu, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(tmp_art_vis.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
    }
  } else {
    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(mu);
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(mu);
    }
  }
  if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
    FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
    tmpMatrix->set_mm_factor(vis_mm_factor.dat);
    tmpMatrix->calc_mat_partial();
  } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
    FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
    tmpMatrix->set_mm_factor(vis_mm_factor.dat);
    tmpMatrix->calc_mat_partial();
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

    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(tmp_art_vis.dat);
    }
  } else {
    if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
      FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(mu);
    } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
      FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
      tmpMatrix->set_factor(mu);
    }
  }
  if(viscositySolverType == LinearSolver::PETSC_JACOBI) {
    FactorMMPoissonMatrixFreeDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeDiag3D*>(viscosityMatrix);
    tmpMatrix->set_mm_factor(vis_mm_factor.dat);
    tmpMatrix->calc_mat_partial();
  } else if(viscositySolverType == LinearSolver::PETSC_BLOCK_JACOBI) {
    FactorMMPoissonMatrixFreeBlockDiag3D *tmpMatrix = dynamic_cast<FactorMMPoissonMatrixFreeBlockDiag3D*>(viscosityMatrix);
    tmpMatrix->set_mm_factor(vis_mm_factor.dat);
    tmpMatrix->calc_mat_partial();
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

  if(values_to_save.count("surface_tension") != 0 && surface_tension) {
    op_fetch_data_hdf5_file(st[(currentInd + 1) % 2][0], filename.c_str());
    op_fetch_data_hdf5_file(st[(currentInd + 1) % 2][1], filename.c_str());
    op_fetch_data_hdf5_file(st[(currentInd + 1) % 2][2], filename.c_str());
  }

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
