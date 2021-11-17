#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

using namespace std;

Poisson::Poisson(INSData *insData, DGMesh *m) {
  data = insData;
  mesh = m;
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
