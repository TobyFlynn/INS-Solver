#include "matrices/3d/factor_mm_poisson_semi_matrix_free_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

FactorMMPoissonSemiMatrixFree3D::FactorMMPoissonSemiMatrixFree3D(DGMesh3D *m) : FactorPoissonSemiMatrixFree3D(m) {

}

void FactorMMPoissonSemiMatrixFree3D::set_mm_factor(op_dat f) {
  mm_factor = f;
}

void FactorMMPoissonSemiMatrixFree3D::calc_mat_partial() {
  timer->startTimer("FactorMMPoissonSemiMatrixFree3D - calc_mat_partial");
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  calc_mm();
  petscMatResetRequired = true;
  timer->endTimer("FactorMMPoissonSemiMatrixFree3D - calc_mat_partial");
}

void FactorMMPoissonSemiMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("FactorMMPoissonSemiMatrixFree3D - Mult");
  FactorPoissonSemiMatrixFree3D::mult(in, out);

  op_par_loop(fpmf_3d_mult_mm, "fpmf_3d_mult_mm", _mesh->cells,
              op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  timer->endTimer("FactorMMPoissonSemiMatrixFree3D - Mult");
  return;
}

void FactorMMPoissonSemiMatrixFree3D::calc_mm() {
  timer->startTimer("FactorMMPoissonSemiMatrixFree3D - calc_mm");
  op_par_loop(factor_poisson_matrix_3d_mm, "factor_poisson_matrix_3d_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("FactorMMPoissonSemiMatrixFree3D - calc_mm");
}
