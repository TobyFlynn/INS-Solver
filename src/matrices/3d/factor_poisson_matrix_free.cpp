#include "matrices/3d/factor_poisson_matrix_free_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

FactorPoissonMatrixFree3D::FactorPoissonMatrixFree3D(DGMesh3D *m) : FactorPoissonSemiMatrixFree3D(m, false) {

}

void FactorPoissonMatrixFree3D::calc_mat() {
  throw std::runtime_error("Cannot calculate explicit matrix of FactorPoissonMatrixFree3D");
}

void FactorPoissonMatrixFree3D::apply_bc(op_dat rhs, op_dat bc) {
  if(mesh->bface2cells) {
    throw std::runtime_error("TODO implement apply_bc of FactorPoissonMatrixFree3D");
  }
}

void FactorPoissonMatrixFree3D::multJacobi(op_dat in, op_dat out) {
  throw std::runtime_error("multJacobi of FactorPoissonMatrixFree3D not implemented");
}

bool FactorPoissonMatrixFree3D::getPETScMat(Mat** mat) {
  throw std::runtime_error("Getting PETSc matrix of FactorPoissonMatrixFree3D not implemented");
  return false;
}

// Empty overrides
void FactorPoissonMatrixFree3D::calc_op1() {

}

void FactorPoissonMatrixFree3D::calc_op2() {

}

void FactorPoissonMatrixFree3D::calc_opbc() {

}

void FactorPoissonMatrixFree3D::calc_glb_ind() {

}
