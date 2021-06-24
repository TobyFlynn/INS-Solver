#include "pressure_solve.h"

#include "op_seq.h"

PressureSolve::PressureSolve(INSData *d, CubatureData *c, GaussData *g) {
  data = d;
  cData = c;
  gData = g;
}

PressureSolve::~PressureSolve() {

}

void PressureSolve::setDirichletBCs(int *d) {
  dirichlet = d;
}

void PressureSolve::setNeumannBCs(int *n) {
  neumann = n;
}

void PressureSolve::setBCValues(op_dat bc) {
  bc_dat = bc;
}

double PressureSolve::getAverageConvergeIter() {
  double res = (double)numberIter/(double)solveCount;
  numberIter = 0;
  solveCount = 0;
  return res;
}

void PressureSolve::init() {

}

bool PressureSolve::solve(op_dat b_dat, op_dat x_dat) {
  return false;
}
