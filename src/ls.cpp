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
  diracDelta_data = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));

  gInput = op_decl_dat(mesh->cells, DG_G_NP, "double", gInput_data, "gInput");

  s      = op_decl_dat(mesh->cells, DG_NP, "double", s_data, "s");
  step_s = op_decl_dat(mesh->cells, DG_NP, "double", step_s_data, "step");
  nx     = op_decl_dat(mesh->cells, DG_NP, "double", nx_data, "ls-nx");
  ny     = op_decl_dat(mesh->cells, DG_NP, "double", ny_data, "ls-ny");
  curv   = op_decl_dat(mesh->cells, DG_NP, "double", curv_data, "curv");
  diracDelta = op_decl_dat(mesh->cells, DG_NP, "double", diracDelta_data, "diracDelta");
}

LS::~LS() {
  free(gInput_data);

  free(s_data);
  free(step_s_data);
  free(nx_data);
  free(ny_data);
  free(curv_data);
  free(diracDelta_data);
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

  F = data->tmp_dg_np[8];
  G = data->tmp_dg_np[9];

  dsdx = data->tmp_dg_np[8];
  dsdy = data->tmp_dg_np[9];

  exAdvec = data->tmp_dg_g_np[0];
  nFlux   = data->tmp_dg_g_np[1];
  gU      = data->tmp_dg_g_np[2];
  gV      = data->tmp_dg_g_np[3];

  op_par_loop(init_surface, "init_surface", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s,       -1, OP_ID, DG_NP, "double", OP_WRITE));

  // h = std::numeric_limits<double>::max();
  h = 0.0;
  op_par_loop(calc_h, "calc_h", mesh->cells,
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(s, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_gbl(&h, 1, "double", OP_MAX));

  // alpha = 2.0 * h / DG_ORDER;
  // order_width = 2.0 * h;
  // epsilon = h / DG_ORDER;
  alpha = 6.0 * h;
  order_width = 6.0 * h;
  epsilon = h;
  reinit_dt = 1.0 / ((DG_ORDER * DG_ORDER / h) + epsilon * ((DG_ORDER * DG_ORDER*DG_ORDER * DG_ORDER)/(h*h)));
  numSteps = ceil((2.0 * alpha / reinit_dt) * 1.1);

  op_par_loop(ls_update_order, "ls_update_order", mesh->cells,
              op_arg_dat(mesh->order,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&order_width,     1, "double", OP_READ),
              op_arg_dat(s,               -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(data->new_order, -1, OP_ID, 1, "int", OP_WRITE));

  std::vector<op_dat> dats_to_update;
  dats_to_update.push_back(data->Q[0][0]);
  dats_to_update.push_back(data->Q[0][1]);
  dats_to_update.push_back(data->Q[1][0]);
  dats_to_update.push_back(data->Q[1][1]);
  dats_to_update.push_back(data->p);
  dats_to_update.push_back(s);

  mesh->update_order(data->new_order, dats_to_update);

  op_par_loop(init_surface, "init_surface", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s,       -1, OP_ID, DG_NP, "double", OP_WRITE));

  // reinit_ls();
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

  // if(reinit_needed())
  //   reinit_ls();

  // Reset LS in boundary cells (this needs is a work around, needs to be fixed properly later)
  // Also has a data race, but all will be setting to same value so should be okay
  /*
  op_par_loop(ls_bcell_reset, "ls_bcell_reset", mesh->bedges,
              op_arg_dat(mesh->order, 0, mesh->bedge2cells, 1, "int", OP_READ),
              op_arg_dat(mesh->x,     0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y,     0, mesh->bedge2cells, DG_NP, "double", OP_READ),
              op_arg_dat(s,           0, mesh->bedge2cells, DG_NP, "double", OP_WRITE));
  */
  /*
  op_par_loop(ls_update_order, "ls_update_order", mesh->cells,
              op_arg_dat(mesh->order,     -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(&order_width,     1, "double", OP_READ),
              op_arg_dat(s,               -1, OP_ID, DG_NP, "double", OP_READ),
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
  update_values();
}

void LS::advec_step(op_dat input, op_dat output) {
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, input, 0.0, gInput);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, u, 0.0, gU);
  op2_gemv(mesh, false, 1.0, DGConstants::GAUSS_INTERP, v, 0.0, gV);

  op_par_loop(zero_g_np, "zero_g_np", mesh->cells,
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
              op_arg_dat(mesh->gauss->x,    0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
              op_arg_dat(mesh->gauss->y,    0, mesh->bedge2cells, DG_G_NP, "double", OP_READ),
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
  op_par_loop(ls_step, "ls_step", mesh->cells,
              op_arg_gbl(&alpha,     1, "double", OP_READ),
              op_arg_dat(s,         -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(step_s,    -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(diracDelta,-1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->rho, -1, OP_ID, DG_NP, "double", OP_WRITE),
              op_arg_dat(data->mu,  -1, OP_ID, DG_NP, "double", OP_WRITE));

  // Assume | grad s | is approx 1 so this is sufficient for getting normals
  grad(mesh, s, nx, ny);

  div(mesh, nx, ny, curv);
}
