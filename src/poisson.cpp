#include "poisson.h"

#include <iostream>
#include <unistd.h>

#include "op_seq.h"

using namespace std;

Poisson::Poisson(DGMesh *m, INSData *d) {
  mesh = m;
  data = d;
}

Poisson::~Poisson() {

}

void Poisson::setDirichletBCs(int *d) {
  dirichlet[0] = d[0];
  dirichlet[1] = d[1];
  dirichlet[2] = d[2];
}

void Poisson::setNeumannBCs(int *n) {
  neumann[0] = n[0];
  neumann[1] = n[1];
  neumann[2] = n[2];
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
