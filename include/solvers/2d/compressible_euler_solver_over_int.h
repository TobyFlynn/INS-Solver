#ifndef __INS_COMPRESSIBLE_EULER_SOLVER_OVER_INT_2D_H
#define __INS_COMPRESSIBLE_EULER_SOLVER_OVER_INT_2D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"

class CompressibleEulerSolverOverInt2D {
public:
  CompressibleEulerSolverOverInt2D(DGMesh2D *m);

  void init();
  void step();
  void dump_data(const std::string &filename);
  DG_FP l2_vortex_error(DG_FP time);

  DG_FP dt;
private:
  void rhs(op_dat *Q, op_dat *RHSQ);

  DGMesh2D *mesh;
  op_dat Q[4], F[4], G[4], gQ[4], gRHSQ[4];
  op_dat rk_wQ[4], rk_RHSQ[3][4];
};

#endif
