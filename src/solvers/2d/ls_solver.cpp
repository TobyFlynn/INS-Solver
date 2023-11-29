#include "solvers/2d/ls_solver.h"

#include "op_seq.h"

#include <limits>
#include <cmath>

#include "dg_op2_blas.h"
#include "dg_constants/dg_constants.h"
#include "dg_dat_pool.h"
#include "timing.h"

extern DGConstants *constants;
extern DGDatPool *dg_dat_pool;
extern Timing *timer;

int counter;

/**************************************************************************
 * LS Advection Solver class that extends the base Advection Solver class *
 **************************************************************************/
LevelSetAdvectionSolver2D::LevelSetAdvectionSolver2D(DGMesh2D *m) : AdvectionSolver2D(m) {}

void LevelSetAdvectionSolver2D::set_bc_types(op_dat bc) {
  bc_types = bc;
}

void LevelSetAdvectionSolver2D::bc_kernel(op_dat val, op_dat u, op_dat v, op_dat out) {
  op_par_loop(ls_advec_2d_bc, "ls_advec_2d_bc", mesh->bfaces,
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

void LevelSetAdvectionSolver2D::bc_kernel_oi(op_dat val, op_dat u, op_dat v, op_dat flux) {
  op_par_loop(ls_advec_2d_oi_bc, "ls_advec_2d_oi_bc", mesh->bfaces,
              op_arg_dat(bc_types,       -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(flux, 0, mesh->bface2cells, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_INC));
}

/************************
 * Main LS Solver class *
 ************************/
LevelSetSolver2D::LevelSetSolver2D(DGMesh2D *m) {
  mesh = m;
  resuming = false;

  s = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_s");
  dsdx = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_dsdx");
  dsdy = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_dsdy");
  s_sample_x = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, (DG_FP *)NULL, "s_sample_x");
  s_sample_y = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, (DG_FP *)NULL, "s_sample_y");

  advecSolver = new LevelSetAdvectionSolver2D(mesh);
}

LevelSetSolver2D::LevelSetSolver2D(DGMesh2D *m, const std::string &filename) {
  mesh = m;
  resuming = true;

  s = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ls_s");

  dsdx = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_dsdx");
  dsdy = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "ls_dsdy");
  s_sample_x = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, (DG_FP *)NULL, "s_sample_x");
  s_sample_y = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, (DG_FP *)NULL, "s_sample_y");

  advecSolver = new LevelSetAdvectionSolver2D(mesh);
}

LevelSetSolver2D::~LevelSetSolver2D() {
  delete advecSolver;
}

void LevelSetSolver2D::init() {
  if(!resuming) {
    op_par_loop(init_surface_2d, "init_surface_2d", mesh->cells,
                op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(s,       -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  }

  // h = std::numeric_limits<DG_FP>::max();
  h = 0.0;
  op_par_loop(calc_h_ls, "calc_h_ls", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));

  op_printf("LS h: %g\n", h);
  // alpha = 2.0 * h / DG_ORDER;
  // order_width = 2.0 * h;
  // epsilon = h / DG_ORDER;
  alpha = 12.0 * h;
  order_width = 12.0 * h;
  epsilon = h;
  // reinit_width = 20.0 * h;
  reinit_width = 50.0;
  reinit_dt = 1.0 / ((DG_ORDER * DG_ORDER / h) + epsilon * ((DG_ORDER * DG_ORDER*DG_ORDER * DG_ORDER)/(h*h)));
  numSteps = ceil((2.0 * alpha / reinit_dt) * 1.1);

  op_printf("Alpha: %g\t\tReinit Width: %g\n", alpha, reinit_width);

  reinitLS();
}

void LevelSetSolver2D::set_bc_types(op_dat bc) {
  advecSolver->set_bc_types(bc);
}

void LevelSetSolver2D::setVelField(op_dat u1, op_dat v1) {
  u = u1;
  v = v1;
}

void LevelSetSolver2D::step(const DG_FP dt, const int num_steps) {
  timer->startTimer("LevelSetSolver2D - step");
  advecSolver->set_dt(dt);
  for(int i = 0; i < num_steps; i++)
    advecSolver->step(s, u, v);

  counter++;
  if(counter > 14) {
    timer->startTimer("LevelSetSolver2D - reinitLS");
    reinitLS();
    timer->endTimer("LevelSetSolver2D - reinitLS");
    counter = 0;
  }
  timer->endTimer("LevelSetSolver2D - step");
}

bool LevelSetSolver2D::reinitNeeded() {
  DG_FP res = 0.0;
  int count = 0;
  mesh->grad(s, dsdx, dsdy);
  op_par_loop(ls_reinit_check, "ls_reinit_check", mesh->cells,
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(s,     -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dsdx,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(dsdy,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(&res,   1, DG_FP_STR, OP_INC),
              op_arg_gbl(&count, 1, "int", OP_INC));

  res = res / (DG_FP)count;
  // std::cout << "LS residual: " << res << " " << abs(1.0 - res) << std::endl;
  return abs(1.0 - res) > 0.01;
}

void LevelSetSolver2D::getRhoMu(op_dat rho, op_dat mu) {
  timer->startTimer("LevelSetSolver2D - getRhoMu");
  op_par_loop(ls_step, "ls_step", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(s,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(mu,  -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("LevelSetSolver2D - getRhoMu");
}

void LevelSetSolver2D::getRhoVolOI(op_dat rho) {
  timer->startTimer("LevelSetSolver2D - getRhoVolOI");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, s, 0.0, rho);
  op_par_loop(ls_step_rho_vol_oi, "ls_step_rho_vol_oi", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));
  timer->endTimer("LevelSetSolver2D - getRhoVolOI");
}

void LevelSetSolver2D::getMuVolOI(op_dat mu) {
  timer->startTimer("LevelSetSolver2D - getMuVolOI");
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, s, 0.0, mu);
  op_par_loop(ls_step_mu_vol_oi, "ls_step_mu_vol_oi", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(mu, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));
  timer->endTimer("LevelSetSolver2D - getMuVolOI");
}

void LevelSetSolver2D::getRhoSurfOI(op_dat rho) {
  timer->startTimer("LevelSetSolver2D - getRhoSurfOI");
  DGTempDat rho_tmp = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  op_par_loop(ls_step_surf_oi_0, "ls_step_surf_oi_0", mesh->cells,
              op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rho_tmp.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, rho_tmp.dat, 0.0, rho);
  dg_dat_pool->releaseTempDatCells(rho_tmp);
  op_par_loop(ls_step_surf_oi_1, "ls_step_surf_oi_1", mesh->cells,
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(rho, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_RW));
  timer->endTimer("LevelSetSolver2D - getRhoSurfOI");
}

void LevelSetSolver2D::getNormalsCurvature(op_dat nx, op_dat ny, op_dat curv) {
  timer->startTimer("LevelSetSolver2D - getNormalsCurvature");
  // Assume | grad s | is approx 1 so this is sufficient for getting normals
  mesh->grad(s, nx, ny);
  mesh->div(nx, ny, curv);
  timer->endTimer("LevelSetSolver2D - getNormalsCurvature");
}

void LevelSetSolver2D::sampleInterface() {
  timer->startTimer("LevelSetSolver2D - sampleInterface");
  op_par_loop(sample_interface, "sample_interface", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(s,           -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(s_sample_x,  -1, OP_ID, LS_SAMPLE_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(s_sample_y,  -1, OP_ID, LS_SAMPLE_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("LevelSetSolver2D - sampleInterface");
}
