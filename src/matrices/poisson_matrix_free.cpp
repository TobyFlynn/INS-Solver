#include "matrices/poisson_matrix_free.h"

#include <stdexcept>

bool PoissonMatrixFree::getPETScMat(Mat** mat) {
  throw std::runtime_error("Not able to get PETSc Matrices for Matrix Free class\n");
  return false;
}

void PoissonMatrixFree::calc_mat() {
  throw std::runtime_error("calc_mat has not been implemented for Matrix Free class\n");
}

void PoissonMatrixFree::mult(op_dat in, op_dat out) {
  throw std::runtime_error("mult has not been implemented for this Matrix Free class\n");
}

void PoissonMatrixFree::multJacobi(op_dat in, op_dat out) {
  throw std::runtime_error("multJacobi has not been implemented for Matrix Free class\n");
}

void PoissonMatrixFree::setPETScMatrix() {
  throw std::runtime_error("Not able to set PETSc Matrices for Matrix Free class\n");
}

void PoissonMatrixFree::calc_op1() {
  throw std::runtime_error("calc_op1 has not been implemented for Matrix Free class\n");
}

void PoissonMatrixFree::calc_op2() {
  throw std::runtime_error("calc_op2 has not been implemented for Matrix Free class\n");
}

void PoissonMatrixFree::calc_opbc() {
  throw std::runtime_error("calc_opbc has not been implemented for Matrix Free class\n");
}

void PoissonMatrixFree::calc_glb_ind() {
  throw std::runtime_error("calc_glb_ind has not been implemented for Matrix Free class\n");
}
