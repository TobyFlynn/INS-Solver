#ifndef __INS_EULER_H
#define __INS_EULER_H

#include <string>

#include "op_seq.h"
#include "dg_mesh.h"

class Euler {
public:
  Euler(std::string &filename);
  ~Euler();

  void step();
  void dump_data(const std::string &filename);

  double dt;
private:
  void init();
  void rhs(op_dat *Q, op_dat *RHSQ);

  DGMesh *mesh;
  op_dat Q[4], F[4], G[4], gQ[4], gRHSQ[4];
  op_dat rk_wQ[4], rk_RHSQ[3][4];
};

#endif
