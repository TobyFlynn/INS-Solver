#include "matrices/3d/poisson_matrix_free_3d.h"

#include "op_seq.h"

#include "dg_constants/dg_constants.h"
#include "dg_op2_blas.h"

#include "timing.h"

extern DGConstants *constants;
extern Timing *timer;

PoissonMatrixFree3D::PoissonMatrixFree3D(DGMesh3D *m) : PoissonSemiMatrixFree3D(m, false) {

}

void PoissonMatrixFree3D::calc_mat() {
  throw std::runtime_error("Cannot calculate explicit matrix of PoissonMatrixFree3D");
}

void PoissonMatrixFree3D::apply_bc(op_dat rhs, op_dat bc) {
  if(mesh->bface2cells) {
    throw std::runtime_error("TODO implement apply_bc of PoissonMatrixFree3D");
  }
}

void PoissonMatrixFree3D::multJacobi(op_dat in, op_dat out) {
  throw std::runtime_error("multJacobi of PoissonMatrixFree3D not implemented");
}

bool PoissonMatrixFree3D::getPETScMat(Mat** mat) {
  throw std::runtime_error("Getting PETSc matrix of PoissonMatrixFree3D not implemented");
  return false;
}

// Empty overrides
void PoissonMatrixFree3D::calc_op1() {

}

void PoissonMatrixFree3D::calc_op2() {

}

void PoissonMatrixFree3D::calc_opbc() {

}

void PoissonMatrixFree3D::calc_glb_ind() {

}
