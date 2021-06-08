#include "ls.h"

#include "op_seq.h"

#include <limits>

#include "constants.h"
#include "blas_calls.h"
#include "operators.h"

#include "kernels/init_surface.h"
#include "kernels/calc_h.h"

#include "kernels/set_rkQ.h"
#include "kernels/update_Q.h"

#include "kernels/ls_advec_edges.h"
#include "kernels/ls_advec_bedges.h"
#include "kernels/ls_advec_flux.h"
#include "kernels/ls_advec_rhs.h"

#include "kernels/ls_sign.h"
#include "kernels/ls_flux.h"
#include "kernels/ls_bflux.h"
#include "kernels/ls_copy.h"
#include "kernels/ls_rhs.h"
#include "kernels/ls_add_diff.h"

LS::LS(INSData *d, CubatureData *c, GaussData *g) {
  data = d;
  cData = c;
  gData = g;

  s_data = (double *)calloc(15 * data->numCells, sizeof(double));

  rk_data[0] = (double *)calloc(15 * data->numCells, sizeof(double));
  rk_data[1] = (double *)calloc(15 * data->numCells, sizeof(double));
  rk_data[2] = (double *)calloc(15 * data->numCells, sizeof(double));
  rkQ_data   = (double *)calloc(15 * data->numCells, sizeof(double));

  F_data       = (double *)calloc(15 * data->numCells, sizeof(double));
  G_data       = (double *)calloc(15 * data->numCells, sizeof(double));
  dFdr_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  dFds_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  dGdr_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  dGds_data    = (double *)calloc(15 * data->numCells, sizeof(double));
  nFlux_data   = (double *)calloc(15 * data->numCells, sizeof(double));
  exAdvec_data = (double *)calloc(15 * data->numCells, sizeof(double));

  dsdx_data   = (double *)calloc(15 * data->numCells, sizeof(double));
  dsdy_data   = (double *)calloc(15 * data->numCells, sizeof(double));
  sign_data   = (double *)calloc(15 * data->numCells, sizeof(double));
  gS_data     = (double *)calloc(21 * data->numCells, sizeof(double));
  dsldx_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  dsrdx_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  dsldy_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  dsrdy_data  = (double *)calloc(21 * data->numCells, sizeof(double));
  dpldx_data  = (double *)calloc(15 * data->numCells, sizeof(double));
  dprdx_data  = (double *)calloc(15 * data->numCells, sizeof(double));
  dpldy_data  = (double *)calloc(15 * data->numCells, sizeof(double));
  dprdy_data  = (double *)calloc(15 * data->numCells, sizeof(double));

  diff_data    = (double *)calloc(15 * data->numCells, sizeof(double));

  s = op_decl_dat(data->cells, 15, "double", s_data, "s");

  rk[0] = op_decl_dat(data->cells, 15, "double", rk_data[0], "rk0");
  rk[1] = op_decl_dat(data->cells, 15, "double", rk_data[1], "rk1");
  rk[2] = op_decl_dat(data->cells, 15, "double", rk_data[2], "rk2");
  rkQ   = op_decl_dat(data->cells, 15, "double", rkQ_data, "rkQ");

  F       = op_decl_dat(data->cells, 15, "double", F_data, "F");
  G       = op_decl_dat(data->cells, 15, "double", G_data, "G");
  dFdr    = op_decl_dat(data->cells, 15, "double", dFdr_data, "dFdr");
  dFds    = op_decl_dat(data->cells, 15, "double", dFds_data, "dFds");
  dGdr    = op_decl_dat(data->cells, 15, "double", dGdr_data, "dGdr");
  dGds    = op_decl_dat(data->cells, 15, "double", dGds_data, "dGds");
  nFlux   = op_decl_dat(data->cells, 15, "double", nFlux_data, "nFlux");
  exAdvec = op_decl_dat(data->cells, 15, "double", exAdvec_data, "exAdvec");

  dsdx   = op_decl_dat(data->cells, 15, "double", dsdx_data, "dsdx");
  dsdy   = op_decl_dat(data->cells, 15, "double", dsdy_data, "dsdy");
  sign   = op_decl_dat(data->cells, 15, "double", sign_data, "sign");
  gS     = op_decl_dat(data->cells, 21, "double", gS_data, "gS");
  dsldx  = op_decl_dat(data->cells, 21, "double", dsldx_data, "dsldx");
  dsrdx  = op_decl_dat(data->cells, 21, "double", dsrdx_data, "dsrdx");
  dsldy  = op_decl_dat(data->cells, 21, "double", dsldy_data, "dsldy");
  dsrdy  = op_decl_dat(data->cells, 21, "double", dsrdy_data, "dsrdy");
  dpldx  = op_decl_dat(data->cells, 15, "double", dpldx_data, "dpldx");
  dprdx  = op_decl_dat(data->cells, 15, "double", dprdx_data, "dprdx");
  dpldy  = op_decl_dat(data->cells, 15, "double", dpldy_data, "dpldy");
  dprdy  = op_decl_dat(data->cells, 15, "double", dprdy_data, "dprdy");

  diff    = op_decl_dat(data->cells, 15, "double", diff_data, "diff");
}

LS::~LS() {
  free(s_data);

  free(rk_data[0]);
  free(rk_data[1]);
  free(rk_data[2]);
  free(rkQ_data);

  free(F_data);
  free(G_data);
  free(dFdr_data);
  free(dFds_data);
  free(dGdr_data);
  free(dGds_data);
  free(nFlux_data);
  free(exAdvec_data);

  free(dsdx_data);
  free(dsdy_data);
  free(sign_data);
  free(gS_data);
  free(dsldx_data);
  free(dsrdx_data);
  free(dsldy_data);
  free(dsrdy_data);
  free(dpldx_data);
  free(dprdx_data);
  free(dpldy_data);
  free(dprdy_data);

  free(diff_data);
}

void LS::init() {
  h = std::numeric_limits<double>::max();
  op_par_loop(calc_h, "calc_h", data->cells,
              op_arg_dat(data->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&h, 1, "double", OP_MIN));

  op_par_loop(init_surface, "init_surface", data->cells,
              op_arg_dat(data->x, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->y, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(s, -1, OP_ID, 15, "double", OP_WRITE));
}

void LS::setVelField(op_dat u1, op_dat v1) {
  u = u1;
  v = v1;
}

void LS::step(double dt) {
  int x = -1;
  op_par_loop(set_rkQ, "set_rkQ", data->cells,
              op_arg_gbl(&x, 1, "int", OP_READ),
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(s, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rk[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rkQ, -1, OP_ID, 15, "double", OP_RW));

  for(int j = 0; j < 3; j++) {
    advec_step(rkQ, rk[j]);

    if(j != 2) {
      op_par_loop(set_rkQ, "set_rkQ", data->cells,
                  op_arg_gbl(&j, 1, "int", OP_READ),
                  op_arg_gbl(&dt, 1, "double", OP_READ),
                  op_arg_dat(s, -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(rk[0], -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(rk[1], -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(rkQ, -1, OP_ID, 15, "double", OP_RW));
    }
  }

  op_par_loop(update_Q, "update_Q", data->cells,
              op_arg_gbl(&dt, 1, "double", OP_READ),
              op_arg_dat(s, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(rk[0], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(rk[2], -1, OP_ID, 15, "double", OP_READ));
}

void LS::advec_step(op_dat input, op_dat output) {
  // Get neighbouring values of q on internal edges
  op_par_loop(ls_advec_edges, "ls_advec_edges", data->edges,
              op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
              op_arg_dat(input, -2, data->edge2cells, 15, "double", OP_READ),
              op_arg_dat(exAdvec, -2, data->edge2cells, 15, "double", OP_INC));

  // Enforce boundary conditions
  op_par_loop(ls_advec_bedges, "ls_advec_bedges", data->bedges,
              op_arg_dat(data->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(data->x, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(data->y, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(input, 0, data->bedge2cells, 15, "double", OP_READ),
              op_arg_dat(exAdvec, 0, data->bedge2cells, 15, "double", OP_INC));

  op_par_loop(ls_advec_flux, "ls_advec_flux", data->cells,
              op_arg_dat(input, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(v, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(F, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(G, -1, OP_ID, 15, "double", OP_WRITE));

  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DRW), 15, F, 0.0, dFdr);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DSW), 15, F, 0.0, dFds);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DRW), 15, G, 0.0, dGdr);
  op2_gemv(true, 15, 15, 1.0, constants->get_ptr(Constants::DSW), 15, G, 0.0, dGds);

  // Calculate vectors F an G from q for each cell
  op_par_loop(ls_advec_rhs, "ls_advec_rhs", data->cells,
              op_arg_dat(dFdr, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dFds, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGdr, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dGds, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->rx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ry, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->sx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->sy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(input, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(exAdvec, -1, OP_ID, 15, "double", OP_RW),
              op_arg_dat(u, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(v, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->fscale, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->nx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(data->ny, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(nFlux, -1, OP_ID, 15, "double", OP_WRITE),
              op_arg_dat(output, -1, OP_ID, 15, "double", OP_WRITE));

  op2_gemv(true, 15, 15, -1.0, constants->get_ptr(Constants::LIFT), 15, nFlux, 1.0, output);
}

void LS::reinit_ls() {
  double alpha = 2.0 * h / 4.0;
  double epsilon = h / 4.0;
  double dt = 1.0 / ((16.0 / h) + epsilon * ((16.0*16.0)/(h*h)));

  grad(data, s, dsdx, dsdy);
  op_par_loop(ls_sign, "ls_sign", data->cells,
              op_arg_gbl(&alpha, 1, "double", OP_READ),
              op_arg_dat(s, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dsdx, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(dsdy, -1, OP_ID, 15, "double", OP_READ),
              op_arg_dat(sign, -1, OP_ID, 15, "double", OP_WRITE));

  double t = 0.0;
  while(t < alpha * 1.5) {
    int x = -1;
    op_par_loop(set_rkQ, "set_rkQ", data->cells,
                op_arg_gbl(&x, 1, "int", OP_READ),
                op_arg_gbl(&dt, 1, "double", OP_READ),
                op_arg_dat(s, -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[0], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[1], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rkQ, -1, OP_ID, 15, "double", OP_RW));

    for(int j = 0; j < 3; j++) {
      op2_gemv(true, 21, 15, 1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, rkQ, 0.0, gS);

      op_par_loop(ls_flux, "ls_flux", data->edges,
                  op_arg_dat(data->edgeNum, -1, OP_ID, 2, "int", OP_READ),
                  op_arg_dat(data->nodeX, -2, data->edge2cells, 3, "double", OP_READ),
                  op_arg_dat(data->nodeY, -2, data->edge2cells, 3, "double", OP_READ),
                  op_arg_dat(gData->sJ, -2, data->edge2cells, 21, "double", OP_READ),
                  op_arg_dat(gData->nx, -2, data->edge2cells, 21, "double", OP_READ),
                  op_arg_dat(gData->ny, -2, data->edge2cells, 21, "double", OP_READ),
                  op_arg_dat(gS, -2, data->edge2cells, 21, "double", OP_READ),
                  op_arg_dat(dsldx, -2, data->edge2cells, 21, "double", OP_INC),
                  op_arg_dat(dsrdx, -2, data->edge2cells, 21, "double", OP_INC),
                  op_arg_dat(dsldy, -2, data->edge2cells, 21, "double", OP_INC),
                  op_arg_dat(dsrdy, -2, data->edge2cells, 21, "double", OP_INC));

      op_par_loop(ls_bflux, "ls_bflux", data->bedges,
                  op_arg_dat(data->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                  op_arg_dat(gData->sJ, 0, data->bedge2cells, 21, "double", OP_READ),
                  op_arg_dat(gData->nx, 0, data->bedge2cells, 21, "double", OP_READ),
                  op_arg_dat(gData->ny, 0, data->bedge2cells, 21, "double", OP_READ),
                  op_arg_dat(gS, 0, data->bedge2cells, 21, "double", OP_READ),
                  op_arg_dat(dsldx, 0, data->bedge2cells, 21, "double", OP_INC),
                  op_arg_dat(dsrdx, 0, data->bedge2cells, 21, "double", OP_INC),
                  op_arg_dat(dsldy, 0, data->bedge2cells, 21, "double", OP_INC),
                  op_arg_dat(dsrdy, 0, data->bedge2cells, 21, "double", OP_INC));

      cub_grad_weak(data, cData, rkQ, dsdx, dsdy);

      op_par_loop(ls_copy, "ls_copy", data->cells,
                  op_arg_dat(dsdx,  -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(dsdy,  -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(dpldx, -1, OP_ID, 15, "double", OP_WRITE),
                  op_arg_dat(dprdx, -1, OP_ID, 15, "double", OP_WRITE),
                  op_arg_dat(dpldy, -1, OP_ID, 15, "double", OP_WRITE),
                  op_arg_dat(dprdy, -1, OP_ID, 15, "double", OP_WRITE));

      op2_gemv(false, 15, 21, -1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, dsldx, 1.0, dpldx);
      op2_gemv(false, 15, 21, -1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, dsrdx, 1.0, dprdx);
      op2_gemv(false, 15, 21, -1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, dsldy, 1.0, dpldy);
      op2_gemv(false, 15, 21, -1.0, constants->get_ptr(Constants::GAUSS_INTERP), 15, dsrdy, 1.0, dprdy);

      op_par_loop(ls_rhs, "ls_rhs", data->cells,
                  op_arg_dat(sign, -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(dpldx, -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(dprdx, -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(dpldy, -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(dprdy, -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(rk[j], -1, OP_ID, 15, "double", OP_WRITE));

      calc_diff(epsilon);

      op_par_loop(ls_add_diff, "ls_add_diff", data->cells,
                  op_arg_dat(diff, -1, OP_ID, 15, "double", OP_READ),
                  op_arg_dat(rk[j], -1, OP_ID, 15, "double", OP_RW),
                  op_arg_dat(dsldx, -1, OP_ID, 21, "double", OP_WRITE),
                  op_arg_dat(dsrdx, -1, OP_ID, 21, "double", OP_WRITE),
                  op_arg_dat(dsldy, -1, OP_ID, 21, "double", OP_WRITE),
                  op_arg_dat(dsrdy, -1, OP_ID, 21, "double", OP_WRITE));

      if(j != 2) {
        op_par_loop(set_rkQ, "set_rkQ", data->cells,
                    op_arg_gbl(&j, 1, "int", OP_READ),
                    op_arg_gbl(&dt, 1, "double", OP_READ),
                    op_arg_dat(s, -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[0], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rk[1], -1, OP_ID, 15, "double", OP_READ),
                    op_arg_dat(rkQ, -1, OP_ID, 15, "double", OP_RW));
      }
    }
    op_par_loop(update_Q, "update_Q", data->cells,
                op_arg_gbl(&dt, 1, "double", OP_READ),
                op_arg_dat(s, -1, OP_ID, 15, "double", OP_RW),
                op_arg_dat(rk[0], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[1], -1, OP_ID, 15, "double", OP_READ),
                op_arg_dat(rk[2], -1, OP_ID, 15, "double", OP_READ));

    t += dt;
  }
}

void LS::calc_diff(double epsilon) {

}
