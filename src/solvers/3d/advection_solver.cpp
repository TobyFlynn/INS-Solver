#include "solvers/3d/advection_solver.h"

#include "op_seq.h"

#include "dg_op2_blas.h"

#include "timing.h"

extern Timing *timer;

AdvectionSolver3D::AdvectionSolver3D(DGMesh3D *m) {
  mesh = m;

  DG_FP *data_t0 = (DG_FP *)calloc(DG_NP * mesh->cells->size, sizeof(DG_FP));
  f = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, data_t0, "advection_f");
  g = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, data_t0, "advection_g");
  h = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, data_t0, "advection_h");
  rk[0] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, data_t0, "advection_rk0");
  rk[1] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, data_t0, "advection_rk1");
  rk[2] = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, data_t0, "advection_rk2");
  rkQ   = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, data_t0, "advection_rkQ");
  free(data_t0);

  DG_FP *data_t1 = (DG_FP *)calloc(DG_NUM_FACES * DG_NPF * mesh->cells->size, sizeof(DG_FP));
  flux = op_decl_dat(mesh->cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, data_t1, "advection_flux");
  free(data_t1);

  dt = -1.0;
}

void AdvectionSolver3D::step(op_dat val, op_dat u, op_dat v, op_dat w) {
  timer->startTimer("AdvectionSolver3D - step");
  if(dt < 0.0)
    set_dt();

  int x = -1;
  op_par_loop(runge_kutta_0, "runge_kutta_0", mesh->cells,
              op_arg_gbl(&x,     1, "int", OP_READ),
              op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
              op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rkQ,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  for(int j = 0; j < 3; j++) {
    rhs(rkQ, u, v, w, rk[j]);

    if(j != 2) {
      op_par_loop(runge_kutta_0, "runge_kutta_0", mesh->cells,
                  op_arg_gbl(&j,     1, "int", OP_READ),
                  op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
                  op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rk[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rk[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
                  op_arg_dat(rkQ,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
    }
  }

  op_par_loop(runge_kutta_1, "runge_kutta_1", mesh->cells,
              op_arg_gbl(&dt,    1, DG_FP_STR, OP_READ),
              op_arg_dat(val,   -1, OP_ID, DG_NP, DG_FP_STR, OP_RW),
              op_arg_dat(rk[0], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(rk[2], -1, OP_ID, DG_NP, DG_FP_STR, OP_READ));
  timer->endTimer("AdvectionSolver3D - step");
}

void AdvectionSolver3D::rhs(op_dat val, op_dat u, op_dat v, op_dat w, op_dat val_out) {
  op_par_loop(advec_3d_0, "advec_3d_0", mesh->cells,
              op_arg_dat(val, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(u,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(v,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(w,   -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(f,   -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(g,   -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE),
              op_arg_dat(h,   -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));

  mesh->div_weak(f, g, h, val_out);

  op_par_loop(zero_npf_1, "zero_npf_1", mesh->cells,
              op_arg_dat(flux, -1, OP_ID, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_WRITE));

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
              op_arg_dat(flux, -2, mesh->face2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));

  if(mesh->bface2cells) {
    op_par_loop(advec_3d_bflux, "advec_3d_bflux", mesh->bfaces,
                op_arg_dat(mesh->bfaceNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->x, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->y, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->z, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bnz, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(val,  0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(u,    0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(v,    0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(w,    0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(flux, 0, mesh->bface2cells, DG_NUM_FACES * DG_NPF, DG_FP_STR, OP_INC));
  }

  op2_gemv(mesh, false, -1.0, DGConstants::LIFT, flux, 1.0, val_out);
}

void AdvectionSolver3D::set_dt() {
  DG_FP h = 0.0;
  op_par_loop(calc_h_3d, "calc_h_3d", mesh->faces,
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_gbl(&h, 1, DG_FP_STR, OP_MAX));
  h = 1.0 / h;
  dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1));
  op_printf("Advection dt is %g\n", dt);
}

void AdvectionSolver3D::set_dt(const DG_FP t) {
  dt = t;
  op_printf("Advection dt is %g\n", dt);
}

void AdvectionSolver3D::set_bc_types(op_dat bc) {
  bc_types = bc;
}
