#include "linear_solvers/linear_solver.h"

#include <stdexcept>

void LinearSolver::set_matrix(PoissonMatrix *mat) {
  matrix = mat;
}

void LinearSolver::set_bcs(op_dat bcs) {
  bc = bcs;
}

void LinearSolver::set_nullspace(bool ns) {
  nullspace = ns;
}

void LinearSolver::init() {

}

void LinearSolver::set_tol(const DG_FP r_tol, const DG_FP a_tol) {
  throw std::runtime_error("set_tol() not implmented");
}
