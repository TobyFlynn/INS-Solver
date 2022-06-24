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
#include "timing.h"

extern Timing *timer;

void LS::reinit_ls() {
  timer->startTimer("LS - Sample Interface");
  op2_gemv(mesh, false, 1.0, DGConstants::DR, s, 0.0, dsdr);
  op2_gemv(mesh, false, 1.0, DGConstants::DS, s, 0.0, dsds);
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

  timer->endTimer("LS - Sample Interface");

  op_arg op2_args[] = {
    op_arg_dat(s_sample_x, -1, OP_ID, LS_SAMPLE_NP, "double", OP_READ),
    op_arg_dat(s_sample_y, -1, OP_ID, LS_SAMPLE_NP, "double", OP_READ),
    op_arg_dat(mesh->x, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->y, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(s, -1, OP_ID, DG_NP, "double", OP_RW),
    op_arg_dat(s_modal, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
    op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
    op_arg_dat(mesh->rx, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->sx, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->ry, -1, OP_ID, DG_NP, "double", OP_READ),
    op_arg_dat(mesh->sy, -1, OP_ID, DG_NP, "double", OP_READ)
  };
  op_mpi_halo_exchanges(s_sample_x->set, 12, op2_args);

  timer->startTimer("LS - Construct K-D Tree");
  KDTree kdtree((double *)s_sample_x->data, (double *)s_sample_y->data, LS_SAMPLE_NP * mesh->numCells);
  timer->endTimer("LS - Construct K-D Tree");

  const double *mesh_x_coords = (double *)mesh->x->data;
  const double *mesh_y_coords = (double *)mesh->y->data;

  timer->startTimer("LS - Newton Method");
  double *closest_x = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  double *closest_y = (double *)calloc(DG_NP * mesh->numCells, sizeof(double));
  int *cell_ind     = (int *)calloc(DG_NP * mesh->numCells, sizeof(int));

  #pragma omp parallel for
  for(int i = 0; i < DG_NP * mesh->numCells; i++) {
    // Get closest sample point
    KDCoord tmp = kdtree.closest_point(mesh_x_coords[i], mesh_y_coords[i]);
    closest_x[i] = tmp.x;
    closest_y[i] = tmp.y;
    cell_ind[i]  = tmp.cell;
  }

  newton_method(DG_NP * mesh->numCells, (double *)s->data, closest_x, closest_y,
                cell_ind, mesh_x_coords, mesh_y_coords, (double *)s_modal->data,
                (double *)mesh->nodeX->data, (double *)mesh->nodeY->data,
                (double *)mesh->rx->data, (double *)mesh->sx->data,
                (double *)mesh->ry->data, (double *)mesh->sy->data, h);

  free(closest_x);
  free(closest_y);
  free(cell_ind);
  timer->endTimer("LS - Newton Method");

  op_mpi_set_dirtybit(12, op2_args);
}
