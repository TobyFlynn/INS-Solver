#include "matrices/3d/factor_mm_poisson_matrix_free_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

void custom_kernel_fpmf_3d_mult_mm(const int order, char const *name, op_set set,
  op_arg arg0,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5);

FactorMMPoissonMatrixFree3D::FactorMMPoissonMatrixFree3D(DGMesh3D *m) : FactorPoissonMatrixFree3D(m) {

}

void FactorMMPoissonMatrixFree3D::set_mm_factor(op_dat f) {
  mm_factor = f;
}

// Doesn't account for BCs
void FactorMMPoissonMatrixFree3D::mult(op_dat in, op_dat out) {
  timer->startTimer("FactorMMPoissonMatrixFree3D - Mult");
  FactorPoissonMatrixFree3D::mult(in, out);

  timer->startTimer("FactorMMPoissonMatrixFree3D - Mult MM");
  #if defined(OP2_DG_CUDA) && !defined(USE_OP2_KERNELS)
  custom_kernel_fpmf_3d_mult_mm(mesh->order_int, "fpmf_3d_mult_mm", _mesh->cells,
              op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  #else
  op_par_loop(fpmf_3d_mult_mm, "fpmf_3d_mult_mm", _mesh->cells,
              op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, DG_FP_STR, OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(in,  -1, OP_ID, DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  #endif
  timer->endTimer("FactorMMPoissonMatrixFree3D - Mult MM");
  timer->endTimer("FactorMMPoissonMatrixFree3D - Mult");
  return;
}
