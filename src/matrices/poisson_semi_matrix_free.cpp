#include "matrices/poisson_semi_matrix_free.h"

#include "op_seq.h"

#include <stdexcept>

#include "timing.h"

extern Timing *timer;

bool PoissonSemiMatrixFree::getPETScMat(Mat** mat) {
  throw std::runtime_error("Not able to get PETSc Matrices for Semi Matrix Free class\n");
  return false;
}

void PoissonSemiMatrixFree::calc_mat() {
  throw std::runtime_error("calc_mat has not been implemented for this Semi Matrix Free class\n");
}

void PoissonSemiMatrixFree::mult(op_dat in, op_dat out) {
  throw std::runtime_error("mult has not been implemented for this Semi Matrix Free class\n");
}

void PoissonSemiMatrixFree::multJacobi(op_dat in, op_dat out) {
  timer->startTimer("PoissonSemiMatrixFree - multJacobi");
  mult(in, out);

  op_par_loop(poisson_mult_jacobi, "poisson_mult_jacobi", _mesh->cells,
              op_arg_dat(_mesh->order, -1, OP_ID, 1, "int", OP_READ),
              op_arg_dat(op1, -1, OP_ID, DG_NP * DG_NP, DG_FP_STR, OP_READ),
              op_arg_dat(out, -1, OP_ID, DG_NP, DG_FP_STR, OP_RW));
  timer->endTimer("PoissonSemiMatrixFree - multJacobi");
}

void PoissonSemiMatrixFree::setPETScMatrix() {
  throw std::runtime_error("Not able to set PETSc Matrices for Semi Matrix Free class\n");
}
