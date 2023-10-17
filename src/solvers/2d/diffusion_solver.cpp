#include "solvers/2d/diffusion_solver.h"

#include "op_seq.h"

#include <limits>
#include <stdexcept>

#include "dg_op2_blas.h"
#include "dg_dat_pool.h"

#include "timing.h"
#include "config.h"

extern Timing *timer;
extern Config *config;
extern DGDatPool *dg_dat_pool;

DiffusionSolver2D::DiffusionSolver2D(DGMesh2D *m) {
  mesh = m;
  dt = -1.0;

  int tmp_oia = 0;
  config->getInt("solver-options", "over_int_advec", tmp_oia);
  over_int_diff = tmp_oia == 1;
}

void DiffusionSolver2D::step(op_dat val, op_dat vis) {
  timer->startTimer("DiffusionSolver2D - step");
  if(dt < 0.0)
    set_dt();

  DGTempDat rk[3];
  rk[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  rk[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  rk[2] = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat rkQ = dg_dat_pool->requestTempDatCells(DG_NP);

  int x = -1;
  op_par_loop(runge_kutta_0, "runge_kutta_0", mesh->cells,
              op_arg_gbl(&x,     1, "int", OP_READ),
              op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
              op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rkQ.dat,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  for(int j = 0; j < 3; j++) {
    if(over_int_diff)
      rhs_over_int(rkQ.dat, vis, rk[j].dat);
    else
      rhs(rkQ.dat, vis, rk[j].dat);

    if(j != 2) {
      op_par_loop(runge_kutta_0, "runge_kutta_0", mesh->cells,
                  op_arg_gbl(&j,     1, "int", OP_READ),
                  op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
                  op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rk[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rk[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rkQ.dat,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
    }
  }

  op_par_loop(runge_kutta_1, "runge_kutta_1", mesh->cells,
              op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
              op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(rk[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  
  dg_dat_pool->releaseTempDatCells(rkQ);
  dg_dat_pool->releaseTempDatCells(rk[0]);
  dg_dat_pool->releaseTempDatCells(rk[1]);
  dg_dat_pool->releaseTempDatCells(rk[2]);
  timer->endTimer("DiffusionSolver2D - step");
}

void DiffusionSolver2D::rhs(op_dat val, op_dat vis, op_dat val_out) {
  DGTempDat val_x = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat val_y = dg_dat_pool->requestTempDatCells(DG_NP);

  mesh->grad_weak(val, val_x.dat, val_y.dat);

  DGTempDat val_flux_x = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat val_flux_y = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(diff_2d_0, "diff_2d_0", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(val, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val_flux_x.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(val_flux_y.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, val_flux_x.dat, -1.0, val_x.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, val_flux_y.dat, -1.0, val_y.dat);

  dg_dat_pool->releaseTempDatCells(val_flux_x);
  dg_dat_pool->releaseTempDatCells(val_flux_y);

  op_par_loop(diff_2d_1, "diff_2d_1", mesh->cells,
              op_arg_dat(vis, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val_x.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(val_y.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  
  mesh->div_weak(val_x.dat, val_y.dat, val_out);

  DGTempDat val_out_flux = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(diff_2d_2, "diff_2d_2", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(val_x.dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val_y.dat, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vis, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val_out_flux.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));
  
  op2_gemv(mesh, false, 1.0, DGConstants::LIFT, val_out_flux.dat, -1.0, val_out);

  dg_dat_pool->releaseTempDatCells(val_out_flux);
  dg_dat_pool->releaseTempDatCells(val_x);
  dg_dat_pool->releaseTempDatCells(val_y);
}

void DiffusionSolver2D::rhs_over_int(op_dat val, op_dat vis, op_dat val_out) {
  throw std::runtime_error("rhs_over_int not implemented for DiffusionSolver2D");
}

void DiffusionSolver2D::set_dt() {
  DG_FP h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  op_printf("h: %g\n", h);

  dt = 1.0 / ((1.0 * DG_ORDER * DG_ORDER / h) + (1.0 * DG_ORDER * DG_ORDER * DG_ORDER * DG_ORDER / (h * h)));
  
  op_printf("dt: %g\n", h);
}

void DiffusionSolver2D::set_dt(const DG_FP t) {
  dt = t;
}
