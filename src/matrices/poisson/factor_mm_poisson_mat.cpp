#include "matrices/2d/factor_mm_poisson_matrix_2d.h"

#include "op_seq.h"

FactorMMPoissonMatrix2D::FactorMMPoissonMatrix2D(DGMesh2D *m) : FactorPoissonMatrix2D(m) {

}

void FactorMMPoissonMatrix2D::set_mm_factor(op_dat f) {
  mm_factor = f;
}

void FactorMMPoissonMatrix2D::calc_mat(op_dat bc_types) {
  calc_glb_ind();
  calc_op1();
  calc_op2();
  calc_opbc(bc_types);
  calc_mm();
  petscMatResetRequired = true;
}

void FactorMMPoissonMatrix2D::calc_mm() {
  op_par_loop(fact_poisson_mm, "fact_poisson_mm", mesh->cells,
              op_arg_dat(mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(mesh->cubature->mm, -1, OP_ID, DG_NP * DG_NP, "double", OP_READ),
              op_arg_dat(mm_factor, -1, OP_ID, DG_NP, "double", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, "double", OP_INC));
}
