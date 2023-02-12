#ifndef __INS_COMPRESSIBLE_EULER_SOLVER_2D_H
#define __INS_COMPRESSIBLE_EULER_SOLVER_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"

class CompressibleEulerSolver2D {
public:
  CompressibleEulerSolver2D(DGMesh2D *m);

  void init();
  void step();
  void dump_data(const std::string &filename);
  double l2_vortex_error(double time);

  double dt;
private:
  void rhs(op_dat *Q, op_dat *RHSQ);

  DGMesh2D *mesh;
  op_dat Q[4], F[4], G[4], gQ[4], gRHSQ[4];
  op_dat rk_wQ[4], rk_RHSQ[3][4];
};

#endif