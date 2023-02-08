#include "solvers/2d/advection_solver.h"

#include "op_seq.h"

#include <limits>

#include "dg_op2_blas.h"

AdvectionSolver2D::AdvectionSolver2D(DGMesh2D *m) {
  mesh = m;

  double *data_t0 = (double *)calloc(DG_NP * mesh->cells->size, sizeof(double));
  f = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advection_2d_f");
  g = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advection_2d_g");
  rk[0] = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advection_2d_rk0");
  rk[1] = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advection_2d_rk1");
  rk[2] = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advection_2d_rk2");
  rkQ   = op_decl_dat(mesh->cells, DG_NP, "double", data_t0, "advection_2d_rkQ");
  free(data_t0);

  double *data_t1 = (double *)calloc(DG_G_NP * mesh->cells->size, sizeof(double));
  flux = op_decl_dat(mesh->cells, DG_G_NP, "double", data_t1, "advection_2d_flux");
  gVal = op_decl_dat(mesh->cells, DG_G_NP, "double", data_t1, "advection_2d_gVal");
  gU   = op_decl_dat(mesh->cells, DG_G_NP, "double", data_t1, "advection_2d_gU");
  gV   = op_decl_dat(mesh->cells, DG_G_NP, "double", data_t1, "advection_2d_gV");
  free(data_t1);

  dt = -1.0;
}

void AdvectionSolver2D::step(op_dat val, op_dat u, op_dat v) {
  if(dt < 0.0)
    set_dt();

  int x = -1;
  op_par_loop(runge_kutta_0, "runge_kutta_0", mesh->cells,
              op_arg_gbl(&x,     1, "int", OP_READ),
              op_arg_gbl(&dt,    1, "double", OP_READ),
              op_arg_dat(val,   -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rkQ,   -1, OP_ID, DG_NP, "double", OP_RW));

  for(int j = 0; j < 3; j++) {
    rhs(rkQ, u, v, rk[j]);

    if(j != 2) {
      op_par_loop(runge_kutta_0, "runge_kutta_0", mesh->cells,
                  op_arg_gbl(&j,     1, "int", OP_READ),
                  op_arg_gbl(&dt,    1, "double", OP_READ),
                  op_arg_dat(val,   -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(rk[0], -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(rk[1], -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(rkQ,   -1, OP_ID, DG_NP, "double", OP_RW));
    }
  }

  op_par_loop(runge_kutta_1, "runge_kutta_1", mesh->cells,
              op_arg_gbl(&dt,    1, "double", OP_READ),
              op_arg_dat(val,   -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(rk[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk[2], -1, OP_ID, DG_NP, "double", OP_READ));
}

void AdvectionSolver2D::rhs(op_dat val, op_dat u, op_dat v, op_dat val_out) {
  op_par_loop(advec_2d_0, "advec_2d_0", mesh->cells,
              op_arg_dat(val, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u,   -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(v,   -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(f,   -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(g,   -1, OP_ID, DG_NP, "double", OP_WRITE));

  mesh->div_weak(f, g, val_out);

  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(flux, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, val, 0.0, gVal);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, u, 0.0, gU);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, v, 0.0, gV);

  op_par_loop(advec_2d_flux, "advec_2d_flux", mesh->faces,
              op_arg_dat(mesh->order,     -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gVal, -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gU,   -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gV,   -2, mesh->face2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(flux, -2, mesh->face2cells, DG_G_NP, "double", OP_INC));

  // TODO BCs

  op2_gemv(mesh, false, -1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, flux, 1.0, val_out);
}

void AdvectionSolver2D::set_dt() {
  double h = std::numeric_limits<double>::max();
  op_par_loop(calc_min_h_2d, "calc_min_h_2d", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&h, 1, "double", OP_MIN));
  dt = h / ((DG_ORDER + 1) * (DG_ORDER + 1));
}

void AdvectionSolver2D::set_dt(const double t) {
  dt = t;
}
