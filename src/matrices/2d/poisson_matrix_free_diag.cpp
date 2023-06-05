#include "matrices/2d/poisson_matrix_free_diag_2d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

PoissonMatrixFreeDiag2D::PoissonMatrixFreeDiag2D(DGMesh2D *m) : PoissonMatrixFreeMult2D(m) {
  _mesh = m;

  diag = op_decl_dat(mesh->cells, DG_NP, DG_FP_STR, (DG_FP *)NULL, "poisson_matrix_free_diag");
}

void PoissonMatrixFreeDiag2D::set_bc_types(op_dat bc_ty) {
  bc_types = bc_ty;
  mat_free_set_bc_types(bc_ty);
}

void PoissonMatrixFreeDiag2D::apply_bc(op_dat rhs, op_dat bc) {
  mat_free_apply_bc(rhs, bc);
}

void PoissonMatrixFreeDiag2D::calc_mat_partial() {
  timer->startTimer("PoissonMatrixFreeDiag2D - calc_mat_partial");
  calc_op1();
  calc_op2();
  calc_opbc();
  petscMatResetRequired = true;
  timer->endTimer("PoissonMatrixFreeDiag2D - calc_mat_partial");
}

void PoissonMatrixFreeDiag2D::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonMatrixFreeDiag2D - mult");
  mat_free_mult(in, out);
  timer->endTimer("PoissonMatrixFreeDiag2D - mult");
}

void PoissonMatrixFreeDiag2D::calc_op1() {
  timer->startTimer("PoissonMatrixFreeDiag2D - calc_op1");
  op_par_loop(poisson_matrix_2d_op1_diag, "poisson_matrix_2d_op1_diag", mesh->cells,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->geof, -1, OP_ID, 5, DG_FP_STR, OP_READ),
              op_arg_dat(diag, -1, OP_ID, DG_NP, DG_FP_STR, OP_WRITE));
  timer->endTimer("PoissonMatrixFreeDiag2D - calc_op1");
}

void PoissonMatrixFreeDiag2D::calc_op2() {
  timer->startTimer("PoissonMatrixFreeDiag2D - calc_op2");
  op_par_loop(poisson_matrix_2d_op2_diag, "poisson_matrix_2d_op2_diag", mesh->faces,
              op_arg_gbl(&mesh->order_int, 1, "int", OP_READ),
              op_arg_dat(mesh->edgeNum, -1, OP_ID, 2, "int", OP_READ),
              op_arg_dat(mesh->reverse, -1, OP_ID, 1, "bool", OP_READ),
              op_arg_dat(mesh->nx, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->ny, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->fscale, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->sJ, -1, OP_ID, 2, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->geof, -2, mesh->face2cells, 5, DG_FP_STR, OP_READ),
              op_arg_dat(diag, 0, mesh->face2cells, DG_NP, DG_FP_STR, OP_INC),
              op_arg_dat(diag, 1, mesh->face2cells, DG_NP, DG_FP_STR, OP_INC));
  timer->endTimer("PoissonMatrixFreeDiag2D - calc_op2");
}

void PoissonMatrixFreeDiag2D::calc_opbc() {
  timer->startTimer("PoissonMatrixFreeDiag2D - calc_opbc");
  // TODO
  timer->endTimer("PoissonMatrixFreeDiag2D - calc_opbc");
}
