#include "matrices/2d/factor_mm_poisson_matrix_over_int_2d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "timing.h"

extern Timing *timer;
extern DGConstants *constants;

FactorMMPoissonMatrixOverInt2D::FactorMMPoissonMatrixOverInt2D(DGMesh2D *m) : FactorPoissonMatrixOverInt2D(m) {

}

void FactorMMPoissonMatrixOverInt2D::set_mm_factor(op_dat f) {
  mm_factor = f;
}

void FactorMMPoissonMatrixOverInt2D::calc_mat() {
  timer->startTimer("FactorMMPoissonMatrixOverInt2D - calc_mat");
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  calc_mm();
  petscMatResetRequired = true;
  timer->endTimer("FactorMMPoissonMatrixOverInt2D - calc_mat");
}

void FactorMMPoissonMatrixOverInt2D::calc_mm() {
  timer->startTimer("FactorMMPoissonMatrixOverInt2D - calc_mm");
  op_par_loop(factor_poisson_matrix_2d_mm, "factor_poisson_matrix_2d_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("FactorMMPoissonMatrixOverInt2D - calc_mm");
}
