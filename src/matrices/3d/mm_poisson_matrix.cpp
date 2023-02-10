#include "matrices/3d/mm_poisson_matrix_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

extern DGConstants *constants;

MMPoissonMatrix3D::MMPoissonMatrix3D(DGMesh3D *m) : PoissonMatrix3D(m) {
  factor = 0.0;
}

void MMPoissonMatrix3D::calc_mat() {
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  calc_mm();
}

void MMPoissonMatrix3D::set_factor(double f) {
  factor = f;
}

double MMPoissonMatrix3D::get_factor() {
  return factor;
}

void MMPoissonMatrix3D::calc_mm() {
  op_par_loop(poisson_matrix_3d_mm, "poisson_matrix_3d_mm", mesh->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_RW));
}
