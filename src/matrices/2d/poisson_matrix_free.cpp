#include "matrices/2d/poisson_matrix_free_2d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

PoissonMatrixFree2D::PoissonMatrixFree2D(DGMesh2D *m) : PoissonMatrixFreeMult2D(m) {
  _mesh = m;
}

void PoissonMatrixFree2D::set_bc_types(op_dat bc_ty) {
  bc_types = bc_ty;
  mat_free_set_bc_types(bc_ty);
}

void PoissonMatrixFree2D::apply_bc(op_dat rhs, op_dat bc) {
  mat_free_apply_bc(rhs, bc);
}

void PoissonMatrixFree2D::mult(op_dat in, op_dat out) {
  timer->startTimer("PoissonMatrixFree2D - mult");
  mat_free_mult(in, out);
  timer->endTimer("PoissonMatrixFree2D - mult");
}
