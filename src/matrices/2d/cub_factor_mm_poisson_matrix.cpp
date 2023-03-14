#include "matrices/2d/cub_factor_mm_poisson_matrix_2d.h"

#include "op_seq.h"

#include "timing.h"

extern Timing *timer;

CubFactorMMPoissonMatrix2D::CubFactorMMPoissonMatrix2D(DGMesh2D *m) : CubFactorPoissonMatrix2D(m) {

}

void CubFactorMMPoissonMatrix2D::set_mm_factor(op_dat f) {
  mm_factor = f;
}

void CubFactorMMPoissonMatrix2D::calc_mat() {
  timer->startTimer("CubFactorMMPoissonMatrix2D - calc_mat");
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  calc_mm();
  petscMatResetRequired = true;
  timer->endTimer("CubFactorMMPoissonMatrix2D - calc_mat");
}

void CubFactorMMPoissonMatrix2D::calc_mm() {
  timer->startTimer("CubFactorMMPoissonMatrix2D - calc_mm");
  op_par_loop(fact_poisson_mm, "fact_poisson_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_INC));
  timer->endTimer("CubFactorMMPoissonMatrix2D - calc_mm");
}
