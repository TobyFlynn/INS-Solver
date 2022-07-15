#include "ls.h"

#include "op_seq.h"

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
  // op2_gemv(mesh, false, 1.0, DGConstants::DR, s, 0.0, dsdr);
  // op2_gemv(mesh, false, 1.0, DGConstants::DS, s, 0.0, dsds);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, s, 0.0, s_modal);
  // op2_gemv(mesh, false, 1.0, DGConstants::INV_V, dsdr, 0.0, dsdr_modal);
  // op2_gemv(mesh, false, 1.0, DGConstants::INV_V, dsds, 0.0, dsds_modal);

  arma::vec x_, y_, r_, s_;
  DGUtils::setRefXY(4, x_, y_);
  DGUtils::xy2rs(x_, y_, r_, s_);

  op_par_loop(sample_interface, "sample_interface", mesh->cells,
              op_arg_gbl(r_.memptr(), LS_SAMPLE_NP, "double", OP_READ),
              op_arg_gbl(s_.memptr(), LS_SAMPLE_NP, "double", OP_READ),
              op_arg_dat(mesh->nodeX, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(mesh->nodeY, -1, OP_ID, 3, "double", OP_READ),
              op_arg_dat(s,           -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(s_modal,     -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->rx,    -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->sx,    -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->ry,    -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(mesh->sy,    -1, OP_ID, DG_NP, "double", OP_READ),
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
  op_mpi_halo_exchanges_cuda(s_sample_x->set, 12, op2_args);

  const int setSize = s_sample_x->set->size;

  double *s_sample_x_h = (double *)calloc(LS_SAMPLE_NP * setSize, sizeof(double));
  cudaMemcpy(s_sample_x_h, s_sample_x->data_d, LS_SAMPLE_NP * setSize * sizeof(double), cudaMemcpyDeviceToHost);
  double *s_sample_y_h = (double *)calloc(LS_SAMPLE_NP * setSize, sizeof(double));
  cudaMemcpy(s_sample_y_h, s_sample_y->data_d, LS_SAMPLE_NP * setSize * sizeof(double), cudaMemcpyDeviceToHost);

  double *mesh_x_h = (double *)calloc(DG_NP * setSize, sizeof(double));
  cudaMemcpy(mesh_x_h, mesh->x->data_d, DG_NP * setSize * sizeof(double), cudaMemcpyDeviceToHost);
  double *mesh_y_h = (double *)calloc(DG_NP * setSize, sizeof(double));
  cudaMemcpy(mesh_y_h, mesh->y->data_d, DG_NP * setSize * sizeof(double), cudaMemcpyDeviceToHost);

  timer->startTimer("LS - Construct K-D Tree");
  KDTree kdtree((double *)s_sample_x_h, (double *)s_sample_y_h, LS_SAMPLE_NP * mesh->numCells);
  timer->endTimer("LS - Construct K-D Tree");

  const double *mesh_x_coords = (double *)mesh_x_h;
  const double *mesh_y_coords = (double *)mesh_y_h;

  timer->startTimer("LS - Newton Method");
  double *closest_x, *closest_y;
  int *cell_ind;
  cudaMallocManaged(&closest_x, DG_NP * mesh->numCells * sizeof(double));
  cudaMallocManaged(&closest_y, DG_NP * mesh->numCells * sizeof(double));
  cudaMallocManaged(&cell_ind, DG_NP * mesh->numCells * sizeof(int));

  #pragma omp parallel for
  for(int i = 0; i < DG_NP * mesh->numCells; i++) {
    // Get closest sample point
    KDCoord tmp = kdtree.closest_point(mesh_x_coords[i], mesh_y_coords[i]);
    closest_x[i] = tmp.x;
    closest_y[i] = tmp.y;
    cell_ind[i]  = tmp.cell;
  }

  newton_method(DG_NP * mesh->numCells, (double *)s->data_d, closest_x,
                closest_y, cell_ind, (double *)mesh->x->data_d,
                (double *)mesh->y->data_d, (double *)s_modal->data_d,
                (double *)mesh->nodeX->data_d, (double *)mesh->nodeY->data_d,
                (double *)mesh->rx->data_d, (double *)mesh->sx->data_d,
                (double *)mesh->ry->data_d, (double *)mesh->sy->data_d, h);

  cudaFree(closest_x);
  cudaFree(closest_y);
  cudaFree(cell_ind);
  timer->endTimer("LS - Newton Method");

  op_mpi_set_dirtybit_cuda(12, op2_args);

  free(s_sample_x_h);
  free(s_sample_y_h);
  free(mesh_x_h);
  free(mesh_y_h);
}
