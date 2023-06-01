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
