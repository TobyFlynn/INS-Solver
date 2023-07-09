#include "matrices/2d/factor_mm_poisson_matrix_free_diag_2d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

FactorMMPoissonMatrixFreeDiag2D::FactorMMPoissonMatrixFreeDiag2D(DGMesh2D *m) : FactorPoissonMatrixFreeDiag2D(m) {

}

void FactorMMPoissonMatrixFreeDiag2D::set_mm_factor(op_dat f) {
  mm_factor = f;
}

void FactorMMPoissonMatrixFreeDiag2D::calc_mat_partial() {
  timer->startTimer("FactorMMPoissonMatrixFreeDiag2D - calc_mat_partial");
  calc_op1();
  calc_op2();
  calc_opbc();
  calc_mm();
  petscMatResetRequired = true;
  timer->endTimer("FactorMMPoissonMatrixFreeDiag2D - calc_mat_partial");
}

void FactorMMPoissonMatrixFreeDiag2D::mult(op_dat in, op_dat out) {
  timer->startTimer("FactorMMPoissonMatrixFreeDiag2D - Mult");
  FactorPoissonMatrixFreeDiag2D::mult(in, out);

  timer->startTimer("FactorMMPoissonMatrixFreeDiag2D - Mult MM");
  op_par_loop(fpmf_2d_mult_mm_geof, "fpmf_2d_mult_mm_geof", _mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("FactorMMPoissonMatrixFreeDiag2D - Mult MM");
  timer->endTimer("FactorMMPoissonMatrixFreeDiag2D - Mult");
}

void FactorMMPoissonMatrixFreeDiag2D::calc_mm() {
  timer->startTimer("FactorMMPoissonMatrixFreeDiag2D - calc_mm");
  op_par_loop(factor_poisson_matrix_2d_mm_diag, "factor_poisson_matrix_2d_mm_diag", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(diag, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("FactorMMPoissonMatrixFreeDiag2D - calc_mm");
}
