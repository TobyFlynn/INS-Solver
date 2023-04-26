#include "matrices/3d/factor_poisson_matrix_free_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

FactorPoissonMatrixFree3D::FactorPoissonMatrixFree3D(DGMesh3D *m) : FactorPoissonMatrixFreeMult3D(m) {

}

void FactorPoissonMatrixFree3D::set_bc_types(op_dat bc_ty) {
  bc_types = bc_ty;
  mat_free_set_bc_types(bc_ty);
}

void FactorPoissonMatrixFree3D::set_factor(op_dat f) {
  factor = f;
  mat_free_set_factor(f);
}

void FactorPoissonMatrixFree3D::apply_bc(op_dat rhs, op_dat bc) {
  mat_free_apply_bc(rhs, bc);
}

void FactorPoissonMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("FactorPoissonMatrixFree3D - mult");
  mat_free_mult(in, out);
  timer->endTimer("FactorPoissonMatrixFree3D - mult");
}
