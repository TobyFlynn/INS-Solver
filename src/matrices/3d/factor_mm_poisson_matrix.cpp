#include "matrices/3d/factor_mm_poisson_matrix_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"

extern DGConstants *constants;

FactorMMPoissonMatrix3D::FactorMMPoissonMatrix3D(DGMesh3D *m) : FactorPoissonMatrix3D(m) {

}

void FactorMMPoissonMatrix3D::set_mm_factor(op_dat f) {
  mm_factor = f;
}

void FactorMMPoissonMatrix3D::calc_mat() {
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc();
  calc_mm();
  petscMatResetRequired = true;
}

void FactorMMPoissonMatrix3D::calc_mm() {
  op_par_loop(factor_poisson_matrix_3d_mm, "factor_poisson_matrix_3d_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_gbl(constants->get_mat_ptr(DGConstants::MASS), DG_ORDER * DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mesh->J, -1, OP_ID, 1, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_RW));
}
