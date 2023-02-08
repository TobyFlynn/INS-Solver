#include "matrices/2d/mm_poisson_matrix_2d.h"

#include "op_seq.h"

MMPoissonMatrix2D::MMPoissonMatrix2D(DGMesh2D *m) : PoissonMatrix2D(m) {
  factor = 0.0;
}

void MMPoissonMatrix2D::calc_mat() {
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc(bc_types);
  calc_mm();
  petscMatResetRequired = true;
}

void MMPoissonMatrix2D::set_factor(double f) {
  factor = f;
}

double MMPoissonMatrix2D::get_factor() {
  return factor;
}

void MMPoissonMatrix2D::calc_mm() {
  op_par_loop(poisson_mm, "poisson_mm", mesh->cells,
              op_arg_gbl(&factor, 1, "double", OP_READ),
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_INC));
}
