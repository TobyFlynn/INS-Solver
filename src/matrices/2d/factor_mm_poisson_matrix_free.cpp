#include "matrices/2d/factor_mm_poisson_matrix_free_2d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

FactorMMPoissonMatrixFree2D::FactorMMPoissonMatrixFree2D(DGMesh2D *m) : FactorPoissonMatrixFree2D(m) {

}

void FactorMMPoissonMatrixFree2D::set_mm_factor(op_dat f) {
  mm_factor = f;
}

// Doesn't account for BCs
void FactorMMPoissonMatrixFree2D::mult(op_dat in, op_dat out) {
  timer->startTimer("FactorMMPoissonMatrixFree2D - Mult");
  FactorPoissonMatrixFree2D::mult(in, out);

  op_par_loop(fpmf_2d_mult_mm, "fpmf_2d_mult_mm", _mesh->cells,
              op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));

  timer->endTimer("FactorMMPoissonMatrixFree2D - Mult");
}
