#ifndef __INS_VISCOSITY_SOLVE_H
#define __INS_VISCOSITY_SOLVE_H

#include "ins_data.h"
#include "op_seq.h"

class ViscositySolve {
public:
  ViscositySolve(INSData *d, CubatureData *c, GaussData *g);
  ~ViscositySolve();

  void init();
  bool solve(op_dat b_dat, op_dat x_dat);

  void setDirichletBCs(int *d);
  void setNeumannBCs(int *n);
  void setBCValues(op_dat bc);
  double getAverageConvergeIter();

private:
  INSData *data;
  CubatureData *cData;
  GaussData *gData;

  int *dirichlet;
  int *neumann;

  int numberIter = 0;
  int solveCount = 0;

  op_dat bc_dat;
};

#endif
