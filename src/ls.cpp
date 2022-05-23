#include "ls.h"

#include "op_seq.h"

#include <limits>
#include <cmath>
#include <vector>

#include <iostream>
#include <fstream>

#include "dg_constants.h"
#include "dg_utils.h"
#include "dg_blas_calls.h"
#include "dg_op2_blas.h"
#include "dg_operators.h"
#include "dg_compiler_defs.h"

#include "kd_tree.h"
#include "utils.h"

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

  s_sample_x_data = (double *)calloc(LS_SAMPLE_NP * mesh->numCells, sizeof(double));
  s_sample_y_data = (double *)calloc(LS_SAMPLE_NP * mesh->numCells, sizeof(double));

  gInput = op_decl_dat(mesh->cells, DG_G_NP, "double", gInput_data, "gInput");

  s      = op_decl_dat(mesh->cells, DG_NP, "double", s_data, "s");
  step_s = op_decl_dat(mesh->cells, DG_NP, "double", step_s_data, "step");
  nx     = op_decl_dat(mesh->cells, DG_NP, "double", nx_data, "ls-nx");
  ny     = op_decl_dat(mesh->cells, DG_NP, "double", ny_data, "ls-ny");
  curv   = op_decl_dat(mesh->cells, DG_NP, "double", curv_data, "curv");
  diracDelta = op_decl_dat(mesh->cells, DG_NP, "double", diracDelta_data, "diracDelta");

  s_sample_x = op_decl_dat(mesh->cells, LS_SAMPLE_NP, "double", s_sample_x_data, "s_sample_x");
  s_sample_y = op_decl_dat(mesh->cells, LS_SAMPLE_NP, "double", s_sample_y_data, "s_sample_y");
}

LS::~LS() {
  free(gInput_data);

  free(s_data);
  free(step_s_data);
  free(nx_data);
  free(ny_data);
  free(curv_data);
  free(diracDelta_data);

  free(s_sample_x_data);
  free(s_sample_y_data);
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

  s_modal     = data->tmp_dg_np[0];
  dsdr_modal  = data->tmp_dg_np[1];
  dsds_modal  = data->tmp_dg_np[2];
  dsdr        = data->tmp_dg_np[3];
  dsds        = data->tmp_dg_np[4];
  dsdr2_modal = data->tmp_dg_np[5];
  dsdrs_modal = data->tmp_dg_np[6];
  dsds2_modal = data->tmp_dg_np[7];

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

  op_printf("h: %g\n", h);
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

  reinit_ls();
  // update_values();
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
  op2_gemv(mesh, false, 1.0, DGConstants::DR, s, 0.0, dsdr);
  op2_gemv(mesh, false, 1.0, DGConstants::DS, s, 0.0, dsds);
  op2_gemv(mesh, false, 1.0, DGConstants::DR, dsdr, 0.0, s_modal);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, s_modal, 0.0, dsdr2_modal);
  op2_gemv(mesh, false, 1.0, DGConstants::DS, dsdr, 0.0, s_modal);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, s_modal, 0.0, dsdrs_modal);
  op2_gemv(mesh, false, 1.0, DGConstants::DS, dsds, 0.0, s_modal);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, s_modal, 0.0, dsds2_modal);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, s, 0.0, s_modal);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, dsdr, 0.0, dsdr_modal);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, dsds, 0.0, dsds_modal);

  arma::vec x_, y_, r_, s_;
  DGUtils::setRefXY(6, x_, y_);
  DGUtils::xy2rs(x_, y_, r_, s_);

  op_par_loop(sample_interface, "sample_interface", mesh->cells,
              op_arg_gbl(r_.memptr(), LS_SAMPLE_NP, "double", OP_READ),
              op_arg_gbl(s_.memptr(), LS_SAMPLE_NP, "double", OP_READ),
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(s,           -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s_modal,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dsdr_modal,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(dsds_modal,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s_sample_x,  -1, OP_ID, LS_SAMPLE_NP, "double", OP_WRITE),
              op_arg_dat(s_sample_y,  -1, OP_ID, LS_SAMPLE_NP, "double", OP_WRITE));

  op_arg op2_args[] = {
    op_arg_dat(s_sample_x, -1, OP_ID, LS_SAMPLE_NP, "double", OP_READ),
    op_arg_dat(s_sample_y, -1, OP_ID, LS_SAMPLE_NP, "double", OP_READ),
    op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(s, -1, OP_ID, DG_NP, "double", OP_RW),
    op_arg_dat(s_modal, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(dsdr_modal,  -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(dsds_modal,  -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(dsdr2_modal, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(dsdrs_modal, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(dsds2_modal, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
    op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
    op_arg_dat(mesh->rx, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->sx, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->ry, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->sy, -1, OP_ID, DG_NP, "double", OP_READ)
  };
  op_mpi_halo_exchanges(s_sample_x->set, 17, op2_args);

  std::string fileName = "interface_sample_points.txt";
  std::ofstream out_file(fileName.c_str());
  out_file << "X,Y,Z" << std::endl;

  const double *x_coords = (double *)s_sample_x->data;
  const double *y_coords = (double *)s_sample_y->data;
  for(int i = 0; i < LS_SAMPLE_NP * mesh->numCells; i++) {
    if(!isnan(x_coords[i]) && !isnan(y_coords[i]))
      out_file << x_coords[i] << "," << y_coords[i] << ",0.0" << std::endl;
  }

  out_file.close();

  KDTree kdtree((double *)s_sample_x->data, (double *)s_sample_y->data, LS_SAMPLE_NP * mesh->numCells);

  const double *mesh_x_coords = (double *)mesh->x->data;
  const double *mesh_y_coords = (double *)mesh->y->data;
  double *dists = (double *)s->data;
  for(int i = 0; i < DG_NP * mesh->numCells; i++) {
    // Get closest sample point
    bool negative = dists[i] < 0.0;
    KDCoord tmp = kdtree.closest_point(mesh_x_coords[i], mesh_y_coords[i]);
    // dists[i] = (tmp.x - mesh_x_coords[i]) * (tmp.x - mesh_x_coords[i]) + (tmp.y - mesh_y_coords[i]) * (tmp.y - mesh_y_coords[i]);
    // dists[i] = sqrt(dists[i]);
    // if(negative) dists[i] *= -1.0;

    // Now use Newton method to get real closest point
    double current_x = tmp.x;
    double current_y = tmp.y;
    const double *s_modal_current_cell = ((double *)s_modal->data) + tmp.cell * DG_NP;
    const double *dsdr_modal_current_cell = ((double *)dsdr_modal->data) + tmp.cell * DG_NP;
    const double *dsds_modal_current_cell = ((double *)dsds_modal->data) + tmp.cell * DG_NP;
    const double *dsdr2_modal_current_cell = ((double *)dsdr2_modal->data) + tmp.cell * DG_NP;
    const double *dsdrs_modal_current_cell = ((double *)dsdrs_modal->data) + tmp.cell * DG_NP;
    const double *dsds2_modal_current_cell = ((double *)dsds2_modal->data) + tmp.cell * DG_NP;
    const double *cellX = ((double *)mesh->nodeX->data) + tmp.cell * 3;
    const double *cellY = ((double *)mesh->nodeY->data) + tmp.cell * 3;
    const double rx = *(((double *)mesh->rx->data) + tmp.cell * DG_NP);
    const double sx = *(((double *)mesh->sx->data) + tmp.cell * DG_NP);
    const double ry = *(((double *)mesh->ry->data) + tmp.cell * DG_NP);
    const double sy = *(((double *)mesh->sy->data) + tmp.cell * DG_NP);
    newton_method(mesh_x_coords[i], mesh_y_coords[i], current_x, current_y,
                  s_modal_current_cell, dsdr_modal_current_cell,
                  dsds_modal_current_cell, dsdr2_modal_current_cell,
                  dsdrs_modal_current_cell, dsds2_modal_current_cell, cellX, cellY, rx, sx, ry, sy);
    dists[i] = (current_x - mesh_x_coords[i]) * (current_x - mesh_x_coords[i]) + (current_y - mesh_y_coords[i]) * (current_y - mesh_y_coords[i]);
    dists[i] = sqrt(dists[i]);
    if(negative) dists[i] *= -1.0;
  }

  op_mpi_set_dirtybit(17, op2_args);

  op_par_loop(ls_reinit_diff, "ls_reinit_diff", mesh->cells,
              op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s,       -1, OP_ID, DG_NP, "double", OP_RW));
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
