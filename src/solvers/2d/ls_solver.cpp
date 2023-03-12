#include "solvers/2d/ls_solver.h"

#include "op_seq.h"

#include <limits>
#include <cmath>

#include "timing.h"

extern Timing *timer;

int counter;

LevelSetSolver2D::LevelSetSolver2D(DGMesh2D *m) {
  mesh = m;
  resuming = false;

  DG_FP *tmp_np = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *tmp_ls_sample_np = (DG_FP *)calloc(LS_SAMPLE_NP * mesh->cells->size, sizeof(DG_FP));

  s = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "ls_s");
  smoothing = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "ls_smoothing"); 
  dsdx = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "ls_dsdx");
  dsdy = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "ls_dsdy");
  s_sample_x = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, tmp_ls_sample_np, "s_sample_x");
  s_sample_y = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, tmp_ls_sample_np, "s_sample_y");

  free(tmp_ls_sample_np);
  free(tmp_np);

  advecSolver = new AdvectionSolver2D(mesh);
}

LevelSetSolver2D::LevelSetSolver2D(DGMesh2D *m, const std::string &filename) {
  mesh = m;
  resuming = true;

  s = op_decl_dat_hdf5(mesh->cells, DG_NP, DG_FP_STR, filename.c_str(), "ls_s");

  DG_FP *tmp_np = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  DG_FP *tmp_ls_sample_np = (DG_FP *)calloc(LS_SAMPLE_NP * mesh->cells->size, sizeof(DG_FP));

  dsdx = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "ls_dsdx");
  dsdy = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, tmp_np, "ls_dsdy");
  s_sample_x = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, tmp_ls_sample_np, "s_sample_x");
  s_sample_y = op_decl_dat(mesh->cells, LS_SAMPLE_NP, DG_FP_STR, tmp_ls_sample_np, "s_sample_y");

  free(tmp_ls_sample_np);
  free(tmp_np);

  advecSolver = new AdvectionSolver2D(mesh);
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
  alpha = 20.0 * h;
  order_width = 12.0 * h;
  epsilon = h;
  reinit_width = 10.0 * h;
  reinit_dt = 1.0 / ((DG_ORDER * DG_ORDER / h) + epsilon * ((DG_ORDER * DG_ORDER*DG_ORDER * DG_ORDER)/(h*h)));
  numSteps = ceil((2.0 * alpha / reinit_dt) * 1.1);
/*
  op_par_loop(ls_update_order, "ls_update_order", mesh->cells,
              op_arg_dat(mesh->order,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&order_width,     1, DG_FP_STR, OP_READ),
              op_arg_dat(s,               -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(data->Q[0][0]);
  dats_to_update.push_back(data->Q[0][1]);
  dats_to_update.push_back(data->Q[1][0]);
  dats_to_update.push_back(data->Q[1][1]);
  dats_to_update.push_back(data->p);
  dats_to_update.push_back(s);

  mesh->update_order(data->new_order, dats_to_update);
*/
  // reinitLS();
}

void LevelSetSolver2D::setVelField(op_dat u1, op_dat v1) {
  u = u1;
  v = v1;
}

void LevelSetSolver2D::step(DG_FP dt) {
  timer->startTimer("LevelSetSolver2D - step");
  advecSolver->set_dt(dt);
  advecSolver->step(s, u, v);

  counter++;
  if(counter > 49) {
    timer->startTimer("LevelSetSolver2D - reinitLS");
    reinitLS();
    timer->endTimer("LevelSetSolver2D - reinitLS");
    counter = 0;
  }

  /*
  op_par_loop(ls_update_order, "ls_update_order", mesh->cells,
              op_arg_dat(mesh->order,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&order_width,     1, DG_FP_STR, OP_READ),
              op_arg_dat(s,               -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(data->Q[0][0]);
  dats_to_update.push_back(data->Q[0][1]);
  dats_to_update.push_back(data->Q[1][0]);
  dats_to_update.push_back(data->Q[1][1]);
  dats_to_update.push_back(data->p);
  dats_to_update.push_back(s);

  mesh->update_order(data->new_order, dats_to_update);
  */
  timer->endTimer("LevelSetSolver2D - step");
}

bool LevelSetSolver2D::reinitNeeded() {
  DG_FP res = 0.0;
  int count = 0;
  mesh->grad(s, dsdx, dsdy);
  op_par_loop(ls_reinit_check, "ls_reinit_check", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
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

void LevelSetSolver2D::getNormalsCurvature(op_dat nx, op_dat ny, op_dat curv) {
  op_par_loop(ls_smoothing, "ls_smoothing", mesh->cells,
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(smoothing, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->cub_grad_with_central_flux(smoothing, nx, ny);

/*
  timer->startTimer("LevelSetSolver2D - getNormalsCurvature");
  // Assume | grad s | is approx 1 so this is sufficient for getting normals
  mesh->cub_grad_with_central_flux(s, nx, ny);
*/

  op_par_loop(ls_normals, "ls_normals", mesh->cells,
              op_arg_dat(s,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(nx, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(ny, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  mesh->cub_div_with_central_flux(nx, ny, curv);

  op_par_loop(ls_normals_curv, "ls_normals_curv", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(s,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(nx, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(ny, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(curv, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
/*
  timer->endTimer("LevelSetSolver2D - getNormalsCurvature");
*/
}

void LevelSetSolver2D::getDiracDelta(op_dat delta_x, op_dat delta_y) {
  op_par_loop(ls_smoothing, "ls_smoothing", mesh->cells,
              op_arg_gbl(&alpha, 1, DG_FP_STR, OP_READ),
              op_arg_dat(s, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(smoothing, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->cub_grad_with_central_flux(smoothing, delta_x, delta_y);
/* 
  timer->startTimer("LevelSetSolver2D - getDiracDelta");
  op_par_loop(ls_delta, "ls_delta", mesh->cells,
              op_arg_gbl(&alpha,  1, DG_FP_STR, OP_READ),
              op_arg_dat(s,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(delta, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("LevelSetSolver2D - getDiracDelta");
*/
}

void LevelSetSolver2D::sampleInterface() {
  timer->startTimer("LevelSetSolver2D - sampleInterface");
  op_par_loop(sample_interface, "sample_interface", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(s,           -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(s_sample_x,  -1, OP_ID, LS_SAMPLE_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(s_sample_y,  -1, OP_ID, LS_SAMPLE_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("LevelSetSolver2D - sampleInterface");
}
