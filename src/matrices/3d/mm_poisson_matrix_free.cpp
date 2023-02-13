#include "matrices/3d/mm_poisson_matrix_free_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

MMPoissonMatrixFree3D::MMPoissonMatrixFree3D(DGMesh3D *m) : PoissonMatrix3D(m) {
  factor = 0.0;
}

void MMPoissonMatrixFree3D::calc_mat() {
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  calc_mm();
}

void MMPoissonMatrixFree3D::set_factor(double f) {
  factor = f;
}

double MMPoissonMatrixFree3D::get_factor() {
  return factor;
}

void MMPoissonMatrixFree3D::calc_mm() {
  op_par_loop(poisson_matrix_3d_mm, "poisson_matrix_3d_mm", mesh->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_RW));
}

// Doesn't account for BCs
void MMPoissonMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("MF - Mult");
  timer->startTimer("MF - Mult - Cells");
  op_par_loop(poisson_mat_free_mult_cells, "poisson_mat_free_mult_cells", _mesh->cells,
              op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->rx, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(mesh->sx, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(mesh->tx, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(mesh->ry, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(mesh->sy, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(mesh->ty, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(mesh->rz, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(mesh->sz, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(mesh->tz, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(mesh->J,  -1, OP_ID, 1, "double", OP_READ),
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, "double", OP_WRITE));
  timer->endTimer("MF - Mult - Cells");
  timer->startTimer("MF - Mult - Faces");
  op_par_loop(poisson_mat_free_mult_faces, "poisson_mat_free_mult_faces", mesh->faces,
              op_arg_dat(mesh->order, -2, mesh->face2cells, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DR), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DS), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::DT), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F0), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F1), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F2), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MM_F3), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->faceNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->fmaskL,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->fmaskR,  -1, OP_ID, DG_NPF, "int", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, "double", OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, "double", OP_READ),
              op_arg_dat(mesh->nz, -1, OP_ID, 2, "double", OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, "double", OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 2, "double", OP_READ),
              op_arg_dat(mesh->rx, -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(mesh->sx, -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(mesh->tx, -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(mesh->ry, -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(mesh->sy, -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(mesh->ty, -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(mesh->rz, -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(mesh->sz, -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(mesh->tz, -2, mesh->face2cells, 1, "double", OP_READ),
              op_arg_dat(in,  -2, mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(out, -2, mesh->face2cells, DG_NP, "double", OP_INC));
  timer->endTimer("MF - Mult - Faces");
  timer->endTimer("MF - Mult");
}
