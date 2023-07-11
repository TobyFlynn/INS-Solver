#ifndef __INS_COMPRESSIBLE_EULER_SOLVER_2D_H
#define __INS_COMPRESSIBLE_EULER_SOLVER_2D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"

#include <vector>

class CompressibleEulerSolver2D {
public:
  CompressibleEulerSolver2D(DGMesh2D *m);

  void init();
  void step();
  void dump_data(const std::string &filename);
  void save_l2_err_history(const std::string &filename);

  DG_FP dt, time;
private:
  void rhs(op_dat *Q, op_dat *RHSQ);
  void record_l2_err();
  DG_FP l2_vortex_error();

  DGMesh2D *mesh;
  op_dat Q[4], F[4], G[4];
  op_dat rk_wQ[4], rk_RHSQ[3][4];
  int l2_counter;

  std::vector<std::pair<DG_FP,DG_FP>> l2_err_history;
};

#endif
