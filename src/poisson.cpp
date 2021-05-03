#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

using namespace std;

Poisson::Poisson(INSData *nsData, CubatureData *cubData, GaussData *gaussData) {
  data = nsData;
  cData = cubData;
  gData = gaussData;
}

Poisson::~Poisson() {

}

void Poisson::setDirichletBCs(int *d) {
  dirichlet = d;
}

void Poisson::setNeumannBCs(int *n) {
  neumann = n;
}

void Poisson::setBCValues(op_dat bc) {
  bc_dat = bc;
}

double Poisson::getAverageConvergeIter() {
  double res = (double)numberIter/(double)solveCount;
  numberIter = 0;
  solveCount = 0;
  return res;
}
