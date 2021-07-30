#include "ls.h"

#include "op_seq.h"

#include <limits>
#include <cmath>
// #include <iostream>

#include "dg_constants.h"
#include "dg_blas_calls.h"
#include "dg_operators.h"
#include "dg_compiler_defs.h"

LS::LS(DGMesh *m, INSData *d) {
  mesh = m;
  data = d;

  s_data      = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  step_s_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  nx_data     = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  ny_data     = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  curv_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  rk_data[0] = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  rk_data[1] = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  rk_data[2] = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  rkQ_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  F_data       = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  G_data       = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  dFdr_data    = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  dFds_data    = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  dGdr_data    = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  dGds_data    = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  nFlux_data   = (double *)calloc(3 * DG_NPF * mesh->numCells, sizeof(double));
  exAdvec_data = (double *)calloc(3 * DG_NPF * mesh->numCells, sizeof(double));

  dsdx_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  dsdy_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  sign_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  gS_data     = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  dsldx_data  = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  dsrdx_data  = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  dsldy_data  = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  dsrdy_data  = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  dpldx_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  dprdx_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  dpldy_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  dprdy_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  sigmax_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  sigmay_data  = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  sigmaFx_data = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  sigmaFy_data = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  gSigmax_data = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  gSigmay_data = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  diff_data    = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  diffF_data   = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));

  s      = op_decl_dat(mesh->cells, DG_NP, "double", s_data, "s");
  step_s = op_decl_dat(mesh->cells, DG_NP, "double", step_s_data, "step");
  nx     = op_decl_dat(mesh->cells, DG_NP, "double", nx_data, "ls-nx");
  ny     = op_decl_dat(mesh->cells, DG_NP, "double", ny_data, "ls-ny");
  curv   = op_decl_dat(mesh->cells, DG_NP, "double", curv_data, "curv");

  rk[0] = op_decl_dat(mesh->cells, DG_NP, "double", rk_data[0], "rk0");
  rk[1] = op_decl_dat(mesh->cells, DG_NP, "double", rk_data[1], "rk1");
  rk[2] = op_decl_dat(mesh->cells, DG_NP, "double", rk_data[2], "rk2");
  rkQ   = op_decl_dat(mesh->cells, DG_NP, "double", rkQ_data, "rkQ");

  F       = op_decl_dat(mesh->cells, DG_NP, "double", F_data, "F");
  G       = op_decl_dat(mesh->cells, DG_NP, "double", G_data, "G");
  dFdr    = op_decl_dat(mesh->cells, DG_NP, "double", dFdr_data, "dFdr");
  dFds    = op_decl_dat(mesh->cells, DG_NP, "double", dFds_data, "dFds");
  dGdr    = op_decl_dat(mesh->cells, DG_NP, "double", dGdr_data, "dGdr");
  dGds    = op_decl_dat(mesh->cells, DG_NP, "double", dGds_data, "dGds");
  nFlux   = op_decl_dat(mesh->cells, 3 * DG_NPF, "double", nFlux_data, "nFlux");
  exAdvec = op_decl_dat(mesh->cells, 3 * DG_NPF, "double", exAdvec_data, "exAdvec");

  dsdx   = op_decl_dat(mesh->cells, DG_NP, "double", dsdx_data, "dsdx");
  dsdy   = op_decl_dat(mesh->cells, DG_NP, "double", dsdy_data, "dsdy");
  sign   = op_decl_dat(mesh->cells, DG_NP, "double", sign_data, "sign");
  gS     = op_decl_dat(mesh->cells, DG_G_NP, "double", gS_data, "gS");
  dsldx  = op_decl_dat(mesh->cells, DG_G_NP, "double", dsldx_data, "dsldx");
  dsrdx  = op_decl_dat(mesh->cells, DG_G_NP, "double", dsrdx_data, "dsrdx");
  dsldy  = op_decl_dat(mesh->cells, DG_G_NP, "double", dsldy_data, "dsldy");
  dsrdy  = op_decl_dat(mesh->cells, DG_G_NP, "double", dsrdy_data, "dsrdy");
  dpldx  = op_decl_dat(mesh->cells, DG_NP, "double", dpldx_data, "dpldx");
  dprdx  = op_decl_dat(mesh->cells, DG_NP, "double", dprdx_data, "dprdx");
  dpldy  = op_decl_dat(mesh->cells, DG_NP, "double", dpldy_data, "dpldy");
  dprdy  = op_decl_dat(mesh->cells, DG_NP, "double", dprdy_data, "dprdy");

  sigmax  = op_decl_dat(mesh->cells, DG_NP, "double", sigmax_data, "sigmax");
  sigmay  = op_decl_dat(mesh->cells, DG_NP, "double", sigmay_data, "sigmay");
  sigmaFx = op_decl_dat(mesh->cells, DG_G_NP, "double", sigmaFx_data, "sigmaFx");
  sigmaFy = op_decl_dat(mesh->cells, DG_G_NP, "double", sigmaFy_data, "sigmaFy");
  gSigmax = op_decl_dat(mesh->cells, DG_G_NP, "double", gSigmax_data, "gSigmax");
  gSigmay = op_decl_dat(mesh->cells, DG_G_NP, "double", gSigmay_data, "gSigmay");
  diff    = op_decl_dat(mesh->cells, DG_NP, "double", diff_data, "diff");
  diffF   = op_decl_dat(mesh->cells, DG_G_NP, "double", diffF_data, "diffF");
}

LS::~LS() {
  free(s_data);
  free(step_s_data);
  free(nx_data);
  free(ny_data);
  free(curv_data);

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

  free(sigmax_data);
  free(sigmay_data);
  free(sigmaFx_data);
  free(sigmaFy_data);
  free(gSigmax_data);
  free(gSigmay_data);
  free(diff_data);
  free(diffF_data);
}

void LS::init() {
  h = std::numeric_limits<double>::max();
  op_par_loop(calc_h, "calc_h", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&h, 1, "double", OP_MIN));

  // alpha = 16.0 * h / 4.0;
  alpha = 8.0 * h;
  epsilon = h / 4.0;
  reinit_dt = 1.0 / ((16.0 / h) + epsilon * ((16.0*16.0)/(h*h)));
  numSteps = ceil((2.0 * alpha / reinit_dt) * 1.1);

  op_par_loop(init_surface, "init_surface", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s,       -1, OP_ID, DG_NP, "double", OP_WRITE));

  update_values();
}

void LS::setVelField(op_dat u1, op_dat v1) {
  u = u1;
  v = v1;
}

void LS::step(double dt) {
  int x = -1;
  op_par_loop(set_rkQ, "set_rkQ", mesh->cells,
              op_arg_gbl(&x,     1, "int", OP_READ),
              op_arg_gbl(&dt,    1, "double", OP_READ),
              op_arg_dat(s,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rkQ,   -1, OP_ID, DG_NP, "double", OP_RW));

  for(int j = 0; j < 3; j++) {
    advec_step(rkQ, rk[j]);

    if(j != 2) {
      op_par_loop(set_rkQ, "set_rkQ", mesh->cells,
                  op_arg_gbl(&j,     1, "int", OP_READ),
                  op_arg_gbl(&dt,    1, "double", OP_READ),
                  op_arg_dat(s,     -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(rk[0], -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(rk[1], -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(rkQ,   -1, OP_ID, DG_NP, "double", OP_RW));
    }
  }

  op_par_loop(update_Q, "update_Q", mesh->cells,
              op_arg_gbl(&dt,    1, "double", OP_READ),
              op_arg_dat(s,     -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(rk[0], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk[1], -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(rk[2], -1, OP_ID, DG_NP, "double", OP_READ));

  if(reinit_needed())
    reinit_ls();

  update_values();
}

void LS::advec_step(op_dat input, op_dat output) {
  // Get neighbouring values of q on internal edges
  op_par_loop(ls_advec_edges, "ls_advec_edges", mesh->edges,
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(input,         -2, mesh->edge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(exAdvec,       -2, mesh->edge2cells, 3 * DG_NPF, "double", OP_INC));

  // Enforce boundary conditions
  op_par_loop(ls_advec_bedges, "ls_advec_bedges", mesh->bedges,
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->x, 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, 0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(input,   0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(exAdvec, 0, mesh->bedge2cells, 3 * DG_NPF, "double", OP_INC));

  op_par_loop(ls_advec_flux, "ls_advec_flux", mesh->cells,
              op_arg_dat(input, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(v,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(F,     -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(G,     -1, OP_ID, DG_NP, "double", OP_WRITE));

  op2_gemv(true, DG_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::DRW), DG_NP, F, 0.0, dFdr);
  op2_gemv(true, DG_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::DSW), DG_NP, F, 0.0, dFds);
  op2_gemv(true, DG_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::DRW), DG_NP, G, 0.0, dGdr);
  op2_gemv(true, DG_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::DSW), DG_NP, G, 0.0, dGds);

  // Calculate vectors F an G from q for each cell
  op_par_loop(ls_advec_rhs, "ls_advec_rhs", mesh->cells,
              op_arg_dat(dFdr,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dFds,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dGdr,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dGds,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->rx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->ry, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->sx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->sy, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(input,    -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(exAdvec,  -1, OP_ID, 3 * DG_NPF, "double", OP_RW),
              op_arg_dat(u,        -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(v,        -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(mesh->nx,     -1, OP_ID, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(mesh->ny,     -1, OP_ID, 3 * DG_NPF, "double", OP_READ),
              op_arg_dat(nFlux,        -1, OP_ID, 3 * DG_NPF, "double", OP_WRITE),
              op_arg_dat(output,       -1, OP_ID, DG_NP, "double", OP_WRITE));

  op2_gemv(true, DG_NP, 3 * DG_NPF, -1.0, constants->get_ptr(DGConstants::LIFT), 3 * DG_NPF, nFlux, 1.0, output);
}

void LS::reinit_ls() {
  cub_grad(mesh, s, dsdx, dsdy);
  inv_mass(mesh, dsdx);
  inv_mass(mesh, dsdy);
  op_par_loop(ls_sign, "ls_sign", mesh->cells,
              op_arg_gbl(&alpha, 1, "double", OP_READ),
              op_arg_dat(s,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dsdx,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dsdy,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(sign,  -1, OP_ID, DG_NP, "double", OP_WRITE));

  for(int i = 0; i < numSteps; i++) {
    int x = -1;
    op_par_loop(set_rkQ, "set_rkQ", mesh->cells,
                op_arg_gbl(&x, 1, "int", OP_READ),
                op_arg_gbl(&reinit_dt, 1, "double", OP_READ),
                op_arg_dat(s,     -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(rk[0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(rk[1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(rkQ,   -1, OP_ID, DG_NP, "double", OP_RW));

    for(int j = 0; j < 3; j++) {
      op2_gemv(true, DG_G_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, rkQ, 0.0, gS);

      op_par_loop(ls_flux, "ls_flux", mesh->edges,
                  op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
                  op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
                  op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
                  op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
                  op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
                  op_arg_dat(gS,    -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
                  op_arg_dat(dsldx, -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
                  op_arg_dat(dsrdx, -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
                  op_arg_dat(dsldy, -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
                  op_arg_dat(dsrdy, -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

      op_par_loop(ls_bflux, "ls_bflux", mesh->bedges,
                  op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                  op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                  op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                  op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                  op_arg_dat(gS,    0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
                  op_arg_dat(dsldx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
                  op_arg_dat(dsrdx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
                  op_arg_dat(dsldy, 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
                  op_arg_dat(dsrdy, 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

      cub_grad_weak(mesh, rkQ, dsdx, dsdy);

      op_par_loop(ls_copy, "ls_copy", mesh->cells,
                  op_arg_dat(dsdx,  -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(dsdy,  -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(dpldx, -1, OP_ID, DG_NP, "double", OP_WRITE),
                  op_arg_dat(dprdx, -1, OP_ID, DG_NP, "double", OP_WRITE),
                  op_arg_dat(dpldy, -1, OP_ID, DG_NP, "double", OP_WRITE),
                  op_arg_dat(dprdy, -1, OP_ID, DG_NP, "double", OP_WRITE));

      op2_gemv(false, DG_NP, DG_G_NP, -1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, dsldx, 1.0, dpldx);
      op2_gemv(false, DG_NP, DG_G_NP, -1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, dsrdx, 1.0, dprdx);
      op2_gemv(false, DG_NP, DG_G_NP, -1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, dsldy, 1.0, dpldy);
      op2_gemv(false, DG_NP, DG_G_NP, -1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, dsrdy, 1.0, dprdy);

      inv_mass(mesh, dpldx);
      inv_mass(mesh, dprdx);
      inv_mass(mesh, dpldy);
      inv_mass(mesh, dprdy);

      op_par_loop(ls_rhs, "ls_rhs", mesh->cells,
                  op_arg_dat(sign,  -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(dpldx, -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(dprdx, -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(dpldy, -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(dprdy, -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(rk[j], -1, OP_ID, DG_NP, "double", OP_WRITE));

      calc_diff();

      op_par_loop(ls_add_diff, "ls_add_diff", mesh->cells,
                  op_arg_dat(diff,  -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(rk[j], -1, OP_ID, DG_NP, "double", OP_RW),
                  op_arg_dat(dsldx, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
                  op_arg_dat(dsrdx, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
                  op_arg_dat(dsldy, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
                  op_arg_dat(dsrdy, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

      if(j != 2) {
        op_par_loop(set_rkQ, "set_rkQ", mesh->cells,
                    op_arg_gbl(&j, 1, "int", OP_READ),
                    op_arg_gbl(&reinit_dt, 1, "double", OP_READ),
                    op_arg_dat(s,     -1, OP_ID, DG_NP, "double", OP_READ),
                    op_arg_dat(rk[0], -1, OP_ID, DG_NP, "double", OP_READ),
                    op_arg_dat(rk[1], -1, OP_ID, DG_NP, "double", OP_READ),
                    op_arg_dat(rkQ,   -1, OP_ID, DG_NP, "double", OP_RW));
      }
    }
    op_par_loop(update_Q, "update_Q", mesh->cells,
                op_arg_gbl(&reinit_dt, 1, "double", OP_READ),
                op_arg_dat(s,     -1, OP_ID, DG_NP, "double", OP_RW),
                op_arg_dat(rk[0], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(rk[1], -1, OP_ID, DG_NP, "double", OP_READ),
                op_arg_dat(rk[2], -1, OP_ID, DG_NP, "double", OP_READ));
  }
}

void LS::calc_diff() {
  // Calculate sigma
  cub_grad_weak(mesh, rkQ, sigmax, sigmay);

  op_par_loop(sigma_flux, "sigma_flux", mesh->edges,
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gS,      -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(sigmaFx, -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(sigmaFy, -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

  op_par_loop(sigma_bflux, "sigma_bflux", mesh->bedges,
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gS,      0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(sigmaFx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(sigmaFy, 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  op2_gemv(false, DG_NP, DG_G_NP, -1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, sigmaFx, 1.0, sigmax);
  op2_gemv(false, DG_NP, DG_G_NP, -1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, sigmaFy, 1.0, sigmay);

  inv_mass(mesh, sigmax);
  inv_mass(mesh, sigmay);

  op_par_loop(sigma_mult, "sigma_mult", mesh->cells,
              op_arg_gbl(&epsilon, 1, "double", OP_READ),
              op_arg_dat(sigmax,  -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(sigmay,  -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(sigmaFx, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(sigmaFy, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(diffF,   -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Calculate diffusion
  cub_div_weak(mesh, sigmax, sigmay, diff);

  op2_gemv(true, DG_G_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, sigmax, 0.0, gSigmax);
  op2_gemv(true, DG_G_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, sigmay, 0.0, gSigmay);

  op_par_loop(diff_flux, "diff_flux", mesh->edges,
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gSigmax, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gSigmay, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(diffF,   -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

  op_par_loop(diff_bflux, "diff_bflux", mesh->bedges,
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gSigmax, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gSigmax, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(diffF,   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  op2_gemv(false, DG_NP, DG_G_NP, -1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, diffF, 1.0, diff);

  inv_mass(mesh, diff);
}

bool LS::reinit_needed() {
  double res = 0.0;
  int count = 0;
  grad(mesh, s, dsdx, dsdy);
  op_par_loop(ls_reinit_check, "ls_reinit_check", mesh->cells,
              op_arg_gbl(&alpha, 1, "double", OP_READ),
              op_arg_dat(s,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dsdx,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dsdy,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_gbl(&res,   1, "double", OP_INC),
              op_arg_gbl(&count, 1, "int", OP_INC));

  res = res / (double)count;
  // std::cout << "LS residual: " << res << " " << abs(1.0 - res) << std::endl;
  return abs(1.0 - res) > 0.1;
}

void LS::update_values() {
  op_par_loop(ls_step, "ls_step", mesh->cells,
              op_arg_gbl(&alpha, 1, "double", OP_READ),
              op_arg_dat(s,         -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(step_s,    -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->nu,  -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_WRITE));

  op2_gemv(true, DG_G_NP, DG_NP, 1.0, constants->get_ptr(DGConstants::GAUSS_INTERP), DG_NP, data->nu, 0.0, data->gNu);

  // Assume | grad s | is approx 1 so this is sufficient for getting normals
  grad(mesh, s, nx, ny);
  div(mesh, nx, ny, curv);
}
