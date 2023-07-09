#include "matrices/2d/factor_poisson_matrix_free_diag_2d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

FactorPoissonMatrixFreeDiag2D::FactorPoissonMatrixFreeDiag2D(DGMesh2D *m) : FactorPoissonMatrixFreeMult2D(m) {
  _mesh = m;

  diag = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "poisson_matrix_free_diag");
}

void FactorPoissonMatrixFreeDiag2D::set_bc_types(op_dat bc_ty) {
  bc_types = bc_ty;
  mat_free_set_bc_types(bc_ty);
}

void FactorPoissonMatrixFreeDiag2D::set_factor(op_dat f) {
  factor = f;
  mat_free_set_factor(f);
}

void FactorPoissonMatrixFreeDiag2D::apply_bc(op_dat rhs, op_dat bc) {
  mat_free_apply_bc(rhs, bc);
}

void FactorPoissonMatrixFreeDiag2D::mult(op_dat in, op_dat out) {
  timer->startTimer("FactorPoissonMatrixFreeDiag2D - mult");
  mat_free_mult(in, out);
  timer->endTimer("FactorPoissonMatrixFreeDiag2D - mult");
}

void FactorPoissonMatrixFreeDiag2D::calc_mat_partial() {
  timer->startTimer("FactorPoissonMatrixFreeDiag2D - calc_mat_partial");
  check_current_order();
  calc_op1();
  calc_op2();
  calc_opbc();
  petscMatResetRequired = true;
  timer->endTimer("FactorPoissonMatrixFreeDiag2D - calc_mat_partial");
}

void FactorPoissonMatrixFreeDiag2D::calc_op1() {
  timer->startTimer("FactorPoissonMatrixFreeDiag2D - calc_op1");
  op_par_loop(factor_poisson_matrix_2d_op1_diag, "factor_poisson_matrix_2d_op1_diag", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_factor_copy, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(diag, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("FactorPoissonMatrixFreeDiag2D - calc_op1");
}

void FactorPoissonMatrixFreeDiag2D::calc_op2() {
  timer->startTimer("FactorPoissonMatrixFreeDiag2D - calc_op2");
  op_par_loop(factor_poisson_matrix_2d_op2_partial_diag, "factor_poisson_matrix_2d_op2_partial_diag", mesh->faces,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -2, mesh->face2cells, 5, DG_FP_STR, OP_READ),
              op_arg_dat(mat_free_factor_copy, -2, mesh->face2cells, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(diag, 0, mesh->face2cells, DG_NP, DG_FP_STR, OP_INC),
              op_arg_dat(diag, 1, mesh->face2cells, DG_NP, DG_FP_STR, OP_INC));
  timer->endTimer("FactorPoissonMatrixFreeDiag2D - calc_op2");
}

void FactorPoissonMatrixFreeDiag2D::calc_opbc() {
  timer->startTimer("FactorPoissonMatrixFreeDiag2D - calc_opbc");
  if(mesh->bface2cells) {
    op_par_loop(factor_poisson_matrix_2d_bop_diag, "factor_poisson_matrix_2d_bop_diag", mesh->bfaces,
                op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
                op_arg_dat(mesh->bedgeNum, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(bc_types, -1, OP_ID, 1, "int", OP_READ),
                op_arg_dat(mesh->bnx, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bny, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bfscale, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->bsJ, -1, OP_ID, 1, DG_FP_STR, OP_READ),
                op_arg_dat(mesh->geof, 0, mesh->bface2cells, 5, DG_FP_STR, OP_READ),
                op_arg_dat(mat_free_factor_copy, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_READ),
                op_arg_dat(diag, 0, mesh->bface2cells, DG_NP, DG_FP_STR, OP_INC));
  }
  timer->endTimer("FactorPoissonMatrixFreeDiag2D - calc_opbc");
}
