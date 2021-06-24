#include "viscosity_solve.h"

#include "op_seq.h"

ViscositySolve::ViscositySolve(INSData *d, CubatureData *c, GaussData *g) {
  data = d;
  cData = c;
  gData = g;
}

ViscositySolve::~ViscositySolve() {

}

void ViscositySolve::setDirichletBCs(int *d) {
  dirichlet = d;
}

void ViscositySolve::setNeumannBCs(int *n) {
  neumann = n;
}

void ViscositySolve::setBCValues(op_dat bc) {
  bc_dat = bc;
}

double ViscositySolve::getAverageConvergeIter() {
  double res = (double)numberIter/(double)solveCount;
  numberIter = 0;
  solveCount = 0;
  return res;
}

void ViscositySolve::init() {

}

bool ViscositySolve::solve(op_dat b_dat, op_dat x_dat) {
  return false;
}
