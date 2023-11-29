#include "solvers/3d/advection_solver.h"

#include "op_seq.h"

#include "dg_op2_blas.h"
#include "dg_dat_pool.h"

#include "config.h"
#include "timing.h"

extern Config *config;
extern Timing *timer;
extern DGDatPool *dg_dat_pool;

AdvectionSolver3D::AdvectionSolver3D(DGMesh3D *m) {
  mesh = m;

  dt = -1.0;
  int tmp_oia = 0;
  config->getInt("solver-options", "over_int_advec", tmp_oia);
  over_int_advec = tmp_oia == 1;
}

void AdvectionSolver3D::step(op_dat val, op_dat u, op_dat v, op_dat w) {
  timer->startTimer("AdvectionSolver3D - step");
  if(dt < 0.0)
    set_dt();

  DGTempDat tmp_rk[3];
  tmp_rk[0] = dg_dat_pool->requestTempDatCells(DG_NP);
  tmp_rk[1] = dg_dat_pool->requestTempDatCells(DG_NP);
  tmp_rk[2] = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_rkQ = dg_dat_pool->requestTempDatCells(DG_NP);

  int x = -1;
  op_par_loop(runge_kutta_0, "runge_kutta_0", mesh->cells,
              op_arg_gbl(&x,     1, "int", OP_READ),
              op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
              op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_rk[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_rk[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_rkQ.dat,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  for(int j = 0; j < 3; j++) {
    if(over_int_advec)
      rhs_over_int(tmp_rkQ.dat, u, v, w, tmp_rk[j].dat);
    else
      rhs(tmp_rkQ.dat, u, v, w, tmp_rk[j].dat);

    if(j != 2) {
      op_par_loop(runge_kutta_0, "runge_kutta_0", mesh->cells,
                  op_arg_gbl(&j,     1, "int", OP_READ),
                  op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
                  op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_rk[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_rk[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(tmp_rkQ.dat,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
    }
  }

  op_par_loop(runge_kutta_1, "runge_kutta_1", mesh->cells,
              op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
              op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_rk[0].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_rk[1].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_rk[2].dat, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));

  dg_dat_pool->releaseTempDatCells(tmp_rk[0]);
  dg_dat_pool->releaseTempDatCells(tmp_rk[1]);
  dg_dat_pool->releaseTempDatCells(tmp_rk[2]);
  dg_dat_pool->releaseTempDatCells(tmp_rkQ);
  timer->endTimer("AdvectionSolver3D - step");
}

void AdvectionSolver3D::rhs(op_dat val, op_dat u, op_dat v, op_dat w, op_dat val_out) {
  DGTempDat tmp_f = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_g = dg_dat_pool->requestTempDatCells(DG_NP);
  DGTempDat tmp_h = dg_dat_pool->requestTempDatCells(DG_NP);

  op_par_loop(advec_3d_0, "advec_3d_0", mesh->cells,
              op_arg_dat(val, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(w,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_f.dat,   -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_g.dat,   -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(tmp_h.dat,   -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div_weak(tmp_f.dat, tmp_g.dat, tmp_h.dat, val_out);
  // mesh->div(tmp_f.dat, tmp_g.dat, tmp_h.dat, val_out);

  dg_dat_pool->releaseTempDatCells(tmp_f);
  dg_dat_pool->releaseTempDatCells(tmp_g);
  dg_dat_pool->releaseTempDatCells(tmp_h);

  DGTempDat tmp_flux = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_NPF);
  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(tmp_flux.dat, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

  op_par_loop(advec_3d_flux, "advec_3d_flux", mesh->faces,
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(val,  -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,    -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,    -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(w,    -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_flux.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));

  if(mesh->bface2cells)
    bc_kernel(val, u, v, w, tmp_flux.dat);

  op2_gemv(mesh, false, -1.0, DGConstants::LIFT, tmp_flux.dat, 1.0, val_out);
  // op2_gemv(mesh, false, 1.0, DGConstants::LIFT, tmp_flux.dat, -1.0, val_out);

  dg_dat_pool->releaseTempDatCells(tmp_flux);
}

void AdvectionSolver3D::rhs_over_int(op_dat val, op_dat u, op_dat v, op_dat w, op_dat val_out) {
  DGTempDat tmp_f = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  DGTempDat tmp_g = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  DGTempDat tmp_h = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);
  DGTempDat val_interp = dg_dat_pool->requestTempDatCells(DG_CUB_3D_NP);

  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, u, 0.0, tmp_f.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, v, 0.0, tmp_g.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, w, 0.0, tmp_h.dat);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_INTERP, val, 0.0, val_interp.dat);

  op_par_loop(advec_3d_oi_0, "advec_3d_oi_0", mesh->cells,
              op_arg_dat(mesh->geof, -1, OP_ID, 10, DG_FP_STR, OP_READ),
              op_arg_dat(val_interp.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_f.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_g.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW),
              op_arg_dat(tmp_h.dat, -1, OP_ID, DG_CUB_3D_NP, DG_FP_STR, OP_RW));

  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PDR, tmp_f.dat, 0.0, val_out);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PDS, tmp_g.dat, 1.0, val_out);
  op2_gemv(mesh, false, 1.0, DGConstants::CUB3D_PDT, tmp_h.dat, 1.0, val_out);

  dg_dat_pool->releaseTempDatCells(tmp_f);
  dg_dat_pool->releaseTempDatCells(tmp_g);
  dg_dat_pool->releaseTempDatCells(tmp_h);
  dg_dat_pool->releaseTempDatCells(val_interp);

  DGTempDat tmp_flux_cub = dg_dat_pool->requestTempDatCells(DG_NUM_FACES * DG_CUB_SURF_3D_NP);
  op_par_loop(zero_cub_surf_3d, "zero_cub_surf_3d", mesh->cells,
              op_arg_dat(tmp_flux_cub.dat, -1, OP_ID, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_WRITE));

  op_par_loop(advec_3d_oi_1, "advec_3d_oi_1", mesh->faces,
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(u,   -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,   -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(w,   -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(val, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(tmp_flux_cub.dat, -2, mesh->face2cells, DG_NUM_FACES * DG_CUB_SURF_3D_NP, DG_FP_STR, OP_WRITE));

  if(mesh->bface2cells)
    bc_kernel_oi(val, u, v, w, tmp_flux_cub.dat);

  op2_gemv(mesh, false, -1.0, DGConstants::CUBSURF3D_LIFT, tmp_flux_cub.dat, 1.0, val_out);
  dg_dat_pool->releaseTempDatCells(tmp_flux_cub);
}

void AdvectionSolver3D::set_dt() {
  DG_FP h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1));
  // op_printf("Advection dt is %g\n", dt);
}

void AdvectionSolver3D::set_dt(const DG_FP t) {
  dt = t;
  // op_printf("Advection dt is %g\n", dt);
}

void AdvectionSolver3D::set_bc_types(op_dat bc) {
  bc_types = bc;
}
