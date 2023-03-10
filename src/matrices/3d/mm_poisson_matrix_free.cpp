#include "matrices/3d/mm_poisson_matrix_free_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

MMPoissonMatrixFree3D::MMPoissonMatrixFree3D(DGMesh3D *m) : PoissonMatrixFree3D(m) {
  factor = 0.0;
}

void MMPoissonMatrixFree3D::set_factor(DG_FP f) {
  factor = f;
}

DG_FP MMPoissonMatrixFree3D::get_factor() {
  return factor;
}

// Doesn't account for BCs
void MMPoissonMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("MMPoissonMatrixFree3D - Mult");
  PoissonMatrixFree3D::mult(in, out);

  op_par_loop(pmf_3d_mult_mm, "pmf_3d_mult_mm", _mesh->cells,
              op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_gbl(&factor,  1, DG_FP_STR, OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  timer->endTimer("MMPoissonMatrixFree3D - Mult");
  return;
}
