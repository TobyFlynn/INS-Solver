#include "ls.h"

#include "op_seq.h"

#include <limits>
#include <cmath>
#include <vector>
// #include <iostream>

#include "dg_constants.h"
#include "dg_blas_calls.h"
#include "dg_op2_blas.h"
#include "dg_operators.h"
#include "dg_compiler_defs.h"

LS::LS(DGMesh *m, INSData *d) {
  mesh = m;
  data = d;

  gInput_data = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));

  s_data      = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  step_s_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  nx_data     = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  ny_data     = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  curv_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  sign_data   = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  gS_data     = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));

  gSigTmp_data = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));
  diff_data    = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  diffF_data   = (double *)calloc(DG_G_NP * mesh->numCells, sizeof(double));

  modal_data     = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  local_vis_data = (double *)calloc(mesh->numCells, sizeof(double));

  gInput = op_decl_dat(mesh->cells, DG_G_NP, "double", gInput_data, "gInput");

  s      = op_decl_dat(mesh->cells, DG_NP, "double", s_data, "s");
  step_s = op_decl_dat(mesh->cells, DG_NP, "double", step_s_data, "step");
  nx     = op_decl_dat(mesh->cells, DG_NP, "double", nx_data, "ls-nx");
  ny     = op_decl_dat(mesh->cells, DG_NP, "double", ny_data, "ls-ny");
  curv   = op_decl_dat(mesh->cells, DG_NP, "double", curv_data, "curv");

  sign   = op_decl_dat(mesh->cells, DG_NP, "double", sign_data, "sign");
  gS     = op_decl_dat(mesh->cells, DG_G_NP, "double", gS_data, "gS");

  gSigTmp = op_decl_dat(mesh->cells, DG_G_NP, "double", gSigTmp_data, "gSigTmp");
  diff    = op_decl_dat(mesh->cells, DG_NP, "double", diff_data, "diff");
  diffF   = op_decl_dat(mesh->cells, DG_G_NP, "double", diffF_data, "diffF");

  modal     = op_decl_dat(mesh->cells, DG_NP, "double", modal_data, "modal");
  local_vis = op_decl_dat(mesh->cells, 1, "double", local_vis_data, "local_vis");
}

LS::~LS() {
  free(gInput_data);

  free(s_data);
  free(step_s_data);
  free(nx_data);
  free(ny_data);
  free(curv_data);

  free(sign_data);
  free(gS_data);

  free(gSigTmp_data);
  free(diff_data);
  free(diffF_data);

  free(modal_data);
  free(local_vis_data);
}

void LS::init() {
  rk[0] = data->tmp_dg_np[0];
  rk[1] = data->tmp_dg_np[1];
  rk[2] = data->tmp_dg_np[2];
  rkQ   = data->tmp_dg_np[3];

  dFdr = data->tmp_dg_np[4];
  dFds = data->tmp_dg_np[5];
  dGdr = data->tmp_dg_np[6];
  dGds = data->tmp_dg_np[7];

  dpldx = data->tmp_dg_np[4];
  dprdx = data->tmp_dg_np[5];
  dpldy = data->tmp_dg_np[6];
  dprdy = data->tmp_dg_np[7];

  F = data->tmp_dg_np[8];
  G = data->tmp_dg_np[9];

  dsdx = data->tmp_dg_np[8];
  dsdy = data->tmp_dg_np[9];

  sigmax = data->tmp_dg_np[8];
  sigmay = data->tmp_dg_np[9];
  sigTmp = data->tmp_dg_np[4];

  exAdvec = data->tmp_dg_g_np[0];
  nFlux   = data->tmp_dg_g_np[1];
  gU      = data->tmp_dg_g_np[2];
  gV      = data->tmp_dg_g_np[3];

  dsldx = data->tmp_dg_g_np[0];
  dsrdx = data->tmp_dg_g_np[1];
  dsldy = data->tmp_dg_g_np[2];
  dsrdy = data->tmp_dg_g_np[3];

  sigmaFx = data->tmp_dg_g_np[0];
  sigmaFy = data->tmp_dg_g_np[1];
  gSigmax = data->tmp_dg_g_np[2];
  gSigmay = data->tmp_dg_g_np[3];

  h = std::numeric_limits<double>::max();
  op_par_loop(calc_h, "calc_h", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_gbl(&h, 1, "double", OP_MIN));

  // alpha = 8.0 * h;
  alpha = 2.0 * h / DG_ORDER;
  epsilon = h / DG_ORDER;
  reinit_dt = 1.0 / ((DG_ORDER * DG_ORDER / h) + epsilon * ((DG_ORDER * DG_ORDER*DG_ORDER * DG_ORDER)/(h*h)));
  numSteps = ceil((2.0 * alpha / reinit_dt) * 1.1);

  op_par_loop(init_surface, "init_surface", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s,       -1, OP_ID, DG_NP, "double", OP_WRITE));

  reinit_ls();
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
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, input, 0.0, gInput);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, u, 0.0, gU);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, v, 0.0, gV);

  op_par_loop(zero_npf, "zero_g_np", mesh->cells,
              op_arg_dat(exAdvec, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(nFlux,   -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Get neighbouring values of q on internal edges
  op_par_loop(ls_advec_edges, "ls_advec_edges", mesh->edges,
              op_arg_dat(mesh->order,   -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(gInput,        -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(exAdvec,       -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

  // Enforce boundary conditions
  op_par_loop(ls_advec_bedges, "ls_advec_bedges", mesh->bedges,
              op_arg_dat(mesh->order, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedge_type, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum,   -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(gInput,  0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(exAdvec, 0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  op_par_loop(ls_advec_flux, "ls_advec_flux", mesh->cells,
              op_arg_dat(input, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(u,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(v,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(F,     -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(G,     -1, OP_ID, DG_NP, "double", OP_WRITE));

  op2_gemv(mesh, false, 1.0, DGConstants::DRW, F, 0.0, dFdr);
  op2_gemv(mesh, false, 1.0, DGConstants::DSW, F, 0.0, dFds);
  op2_gemv(mesh, false, 1.0, DGConstants::DRW, G, 0.0, dGdr);
  op2_gemv(mesh, false, 1.0, DGConstants::DSW, G, 0.0, dGds);

  // Calculate vectors F an G from q for each cell
  op_par_loop(ls_advec_rhs, "ls_advec_rhs", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(dFdr,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dFds,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dGdr,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dGds,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->rx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->ry, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->sx, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->sy, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(gInput,   -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(exAdvec,  -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(gU,       -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(gV,       -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -1, OP_ID, DG_G_NP, "double", OP_READ),
              op_arg_dat(nFlux,  -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(output, -1, OP_ID, DG_NP, "double", OP_WRITE));

  op2_gemv(mesh, false, -1.0, DGConstants::INV_MASS_GAUSS_INTERP_T, nFlux, 1.0, output);
}

void LS::reinit_ls() {
  // std::cout << "***reinit ls***" << std::endl;
  calc_local_diff_const();

  // cub_grad(mesh, s, dsdx, dsdy);
  grad(mesh, s, dsdx, dsdy);
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
      op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, rkQ, 0.0, gS);

      op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
                  op_arg_dat(dsldx, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
                  op_arg_dat(dsrdx, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

      op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
                  op_arg_dat(dsldy, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
                  op_arg_dat(dsrdy, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

      op_par_loop(ls_flux, "ls_flux", mesh->edges,
                  op_arg_dat(mesh->order,     -2, mesh->edge2cells, 1, "int", OP_READ),
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
                  op_arg_dat(mesh->order,     0, mesh->bedge2cells, 1, "int", OP_READ),
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
      // grad_weak(mesh, rkQ, dsdx, dsdy);

      op_par_loop(ls_copy, "ls_copy", mesh->cells,
                  op_arg_dat(dsdx,  -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(dsdy,  -1, OP_ID, DG_NP, "double", OP_READ),
                  op_arg_dat(dpldx, -1, OP_ID, DG_NP, "double", OP_WRITE),
                  op_arg_dat(dprdx, -1, OP_ID, DG_NP, "double", OP_WRITE),
                  op_arg_dat(dpldy, -1, OP_ID, DG_NP, "double", OP_WRITE),
                  op_arg_dat(dprdy, -1, OP_ID, DG_NP, "double", OP_WRITE));

      op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, dsldx, 1.0, dpldx);
      op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, dsrdx, 1.0, dprdx);
      op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, dsldy, 1.0, dpldy);
      op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, dsrdy, 1.0, dprdy);

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
                  op_arg_dat(rk[j], -1, OP_ID, DG_NP, "double", OP_RW));

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
  op_par_loop(sigma_mult, "sigma_mult", mesh->cells,
              op_arg_dat(local_vis, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(rkQ,       -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(sigTmp,    -1, OP_ID, DG_NP, "double", OP_WRITE));

  // Calculate sigma
  cub_grad_weak(mesh, sigTmp, sigmax, sigmay);
  // grad_weak(mesh, sigTmp, sigmax, sigmay);

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
              op_arg_dat(sigmaFx, -1, OP_ID, DG_G_NP, "double", OP_WRITE),
              op_arg_dat(sigmaFy, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  op_par_loop(sigma_flux, "sigma_flux", mesh->edges,
              op_arg_dat(mesh->order,     -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gS,        -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(local_vis, -2, mesh->edge2cells, 1, "double", OP_READ),
              op_arg_dat(sigmaFx,   -2, mesh->edge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(sigmaFy,   -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

  op_par_loop(sigma_bflux, "sigma_bflux", mesh->bedges,
              op_arg_dat(mesh->order,     0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gS,        0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(local_vis, 0, mesh->bedge2cells, 1, "double", OP_READ),
              op_arg_dat(sigmaFx,   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC),
              op_arg_dat(sigmaFy,   0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, sigmaFx, 1.0, sigmax);
  op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, sigmaFy, 1.0, sigmay);

  inv_mass(mesh, sigmax);
  inv_mass(mesh, sigmay);

  op_par_loop(zero_g_np1, "zero_g_np1", mesh->cells,
              op_arg_dat(diffF, -1, OP_ID, DG_G_NP, "double", OP_WRITE));

  // Calculate diffusion
  op_par_loop(diff_mult, "diff_mult", mesh->cells,
              op_arg_dat(local_vis, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(sigmax,    -1, OP_ID, DG_NP, "double", OP_RW),
              op_arg_dat(sigmay,    -1, OP_ID, DG_NP, "double", OP_RW));

  cub_div_weak(mesh, sigmax, sigmay, diff);
  // div_weak(mesh, sigmax, sigmay, diff);

  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, sigmax, 0.0, gSigmax);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, sigmay, 0.0, gSigmay);

  op_par_loop(diff_flux, "diff_flux", mesh->edges,
              op_arg_dat(mesh->order,     -2, mesh->edge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum,   -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse,   -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->gauss->sJ, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gS,        -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(local_vis, -2, mesh->edge2cells, 1, "double", OP_READ),
              op_arg_dat(gSigmax,   -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gSigmay,   -2, mesh->edge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(diffF,     -2, mesh->edge2cells, DG_G_NP, "double", OP_INC));

  op_par_loop(diff_bflux, "diff_bflux", mesh->bedges,
              op_arg_dat(mesh->order,     0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->gauss->sJ, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->nx, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->ny, 0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gS,        0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(local_vis, 0, mesh->bedge2cells, 1, "double", OP_READ),
              op_arg_dat(gSigmax,   0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(gSigmax,   0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(diffF,     0, mesh->bedge2cells, DG_G_NP, "double", OP_INC));

  op2_gemv(mesh, true, -1.0, DGConstants::GAUSS_INTERP, diffF, 1.0, diff);

  inv_mass(mesh, diff);
}

bool LS::reinit_needed() {
  double res = 0.0;
  int count = 0;
  grad(mesh, s, dsdx, dsdy);
  op_par_loop(ls_reinit_check, "ls_reinit_check", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&alpha, 1, "double", OP_READ),
              op_arg_dat(s,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dsdx,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dsdy,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_gbl(&res,   1, "double", OP_INC),
              op_arg_gbl(&count, 1, "int", OP_INC));

  res = res / (double)count;
  // std::cout << "LS residual: " << res << " " << abs(1.0 - res) << std::endl;
  return abs(1.0 - res) > 0.01;
}

void LS::update_values() {
  op_par_loop(ls_update_order, "ls_update_order", mesh->cells,
              op_arg_dat(mesh->order,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&alpha,           1, "double", OP_READ),
              op_arg_dat(s,               -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(data->Q[0][0]);
  dats_to_update.push_back(data->Q[0][1]);
  dats_to_update.push_back(data->Q[1][0]);
  dats_to_update.push_back(data->Q[1][1]);
  dats_to_update.push_back(s);

  mesh->update_order(data->new_order, dats_to_update);


  op_par_loop(ls_step, "ls_step", mesh->cells,
              op_arg_gbl(&alpha,  1, "double", OP_READ),
              op_arg_dat(s,      -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(step_s, -1, OP_ID, DG_NP, "double", OP_WRITE));

  // Assume | grad s | is approx 1 so this is sufficient for getting normals
  grad(mesh, s, nx, ny);

  div(mesh, nx, ny, curv);
}

void LS::calc_local_diff_const() {
  // Get modal coefficients from nodal representation
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, s, 0.0, modal);

  // Group modal coefficients using quadratic mean (with skyline pessimization)
  op_par_loop(ls_local_vis, "ls_local_vis", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&epsilon,   1, "double", OP_READ),
              op_arg_dat(modal,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(local_vis, -1, OP_ID, 1, "double", OP_WRITE));
}
