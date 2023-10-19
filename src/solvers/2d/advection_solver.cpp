#include "solvers/2d/advection_solver.h"

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

AdvectionSolver2D::AdvectionSolver2D(DGMesh2D *m) {
  mesh = m;
  dt = -1.0;

  int tmp_oia = 0;
  config->getInt("solver-options", "over_int_advec", tmp_oia);
  over_int_advec = tmp_oia == 1;
}

void AdvectionSolver2D::step(op_dat val, op_dat u, op_dat v) {
  timer->startTimer("AdvectionSolver2D - step");
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
    if(over_int_advec)
      rhs_over_int(rkQ.dat, u, v, rk[j].dat);
    else
      rhs(rkQ.dat, u, v, rk[j].dat);

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
  timer->endTimer("AdvectionSolver2D - step");
}

void AdvectionSolver2D::rhs(op_dat val, op_dat u, op_dat v, op_dat val_out) {
  DGTempDat f = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat g = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(advec_2d_0, "advec_2d_0", mesh->cells,
              op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,     -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,     -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(g.dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div_weak(f.dat, g.dat, val_out);

  dg_dat_pool->releaseTempDatCells(f);
  dg_dat_pool->releaseTempDatCells(g);

  DGTempDat flux = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(flux.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(advec_2d_flux, "advec_2d_flux", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(val, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,   -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,   -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(flux.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));

  op_par_loop(advec_2d_bflux, "advec_2d_bflux", mesh->bfaces,
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,   0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(flux.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));

  op2_gemv(mesh, false, -1.0, DGConstants::LIFT, flux.dat, 1.0, val_out);

  dg_dat_pool->releaseTempDatCells(flux);
}

void AdvectionSolver2D::rhs_over_int(op_dat val, op_dat u, op_dat v, op_dat val_out) {
  DGTempDat f = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat g = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);
  DGTempDat val_interp = dg_dat_pool->requestTempDatCells(DG_CUB_2D_NP);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, u, 0.0, f.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, v, 0.0, g.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_INTERP, val, 0.0, val_interp.dat);

  op_par_loop(advec_2d_oi_0, "advec_2d_oi_0", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(val_interp.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(g.dat, -1, OP_ID, DG_CUB_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDR, f.dat, 0.0, val_out);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB2D_PDS, g.dat, 1.0, val_out);

  dg_dat_pool->releaseTempDatCells(f);
  dg_dat_pool->releaseTempDatCells(g);
  dg_dat_pool->releaseTempDatCells(val_interp);

  DGTempDat uM = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat vM = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat valM = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  DGTempDat valP = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);

  op_par_loop(advec_2d_oi_1, "advec_2d_oi_1", mesh->faces,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(u, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(uM.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(vM.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(valM.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE),
              op_arg_dat(valP.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  if(mesh->bface2cells) {
    op_par_loop(advec_2d_oi_2, "advec_2d_oi_2", mesh->bfaces,
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(u, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(val, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(uM.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
                op_arg_dat(vM.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
                op_arg_dat(valM.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW),
                op_arg_dat(valP.dat, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_RW));
  }

  DGTempDat uM_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  DGTempDat vM_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  DGTempDat valM_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);
  DGTempDat valP_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_2D_NP);

  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, uM.dat, 0.0, uM_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, vM.dat, 0.0, vM_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, valM.dat, 0.0, valM_cub.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUBSURF2D_INTERP, valP.dat, 0.0, valP_cub.dat);

  dg_dat_pool->releaseTempDatCells(uM);
  dg_dat_pool->releaseTempDatCells(vM);
  dg_dat_pool->releaseTempDatCells(valM);
  dg_dat_pool->releaseTempDatCells(valP);

  op_par_loop(advec_2d_oi_3, "advec_2d_oi_3", mesh->cells,
              op_arg_dat(mesh->nx_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ_c, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(uM_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(vM_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(valM_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(valP_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_2D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, -1.0, DGConstants::CUBSURF2D_LIFT, valP_cub.dat, 1.0, val_out);

  dg_dat_pool->releaseTempDatCells(uM_cub);
  dg_dat_pool->releaseTempDatCells(vM_cub);
  dg_dat_pool->releaseTempDatCells(valM_cub);
  dg_dat_pool->releaseTempDatCells(valP_cub);
}

void AdvectionSolver2D::set_dt() {
  DG_FP h = std::numeric_limits<DG_FP>::max();
  op_par_loop(calc_min_h_2d, "calc_min_h_2d", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MIN));
  dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1));
}

void AdvectionSolver2D::set_dt(const DG_FP t) {
  dt = t;
}
