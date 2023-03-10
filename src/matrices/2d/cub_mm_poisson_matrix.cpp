#include "matrices/2d/cub_mm_poisson_matrix_2d.h"

#include "op_seq.h"

#include "timing.h"

extern Timing *timer;

CubMMPoissonMatrix2D::CubMMPoissonMatrix2D(DGMesh2D *m) : CubPoissonMatrix2D(m) {
  factor = 0.0;
}

void CubMMPoissonMatrix2D::calc_mat() {
  timer->startTimer("CubMMPoissonMatrix2D - calc_mat");
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  calc_mm();
  petscMatResetRequired = true;
  timer->endTimer("CubMMPoissonMatrix2D - calc_mat");
}

void CubMMPoissonMatrix2D::set_factor(DG_FP f) {
  factor = f;
}

DG_FP CubMMPoissonMatrix2D::get_factor() {
  return factor;
}

void CubMMPoissonMatrix2D::calc_mm() {
  timer->startTimer("CubMMPoissonMatrix2D - calc_mm");
  op_par_loop(poisson_mm, "poisson_mm", mesh->cells,
              op_arg_gbl(&factor, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_INC));
  timer->endTimer("CubMMPoissonMatrix2D - calc_mm");
}
