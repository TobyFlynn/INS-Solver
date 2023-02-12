#include "matrices/3d/mm_poisson_matrix_free_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

extern DGConstants *constants;

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

  op_par_loop(poisson_mult_faces, "poisson_mult_faces", _mesh->faces,
              op_arg_dat(_mesh->order, 0, _mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(in,      0, _mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[0], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out,     0, _mesh->face2cells, DG_NP, "double", OP_INC),
              op_arg_dat(_mesh->order, 1, _mesh->face2cells, 1, "int", OP_READ),
              op_arg_dat(in,      1, _mesh->face2cells, DG_NP, "double", OP_READ),
              op_arg_dat(op2[1], -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(out,     1, _mesh->face2cells, DG_NP, "double", OP_INC));
}