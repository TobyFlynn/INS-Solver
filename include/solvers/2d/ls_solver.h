#ifndef __INS_LS_2D_H
#define __INS_LS_2D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "solvers/2d/advection_solver.h"

#include "dg_mesh/dg_mesh_2d.h"

class LevelSetSolver2D {
public:
  LevelSetSolver2D(DGMesh2D *m);
  ~LevelSetSolver2D();

  void init();

  void setVelField(op_dat u1, op_dat v1);
  void step(DG_FP dt);
  void getRhoMu(op_dat rho, op_dat mu);
  void getNormalsCurvature(op_dat nx, op_dat ny, op_dat curv);
  void getDiracDelta(op_dat delta);

  DGMesh2D *mesh;

  op_dat u, v, s, dsdx, dsdy, s_sample_x, s_sample_y;

  DG_FP alpha, order_width;
private:
  void sampleInterface();
  void reinitLS();
  bool reinitNeeded();

  DG_FP h, epsilon, reinit_dt, reinit_width;
  int numSteps;

  AdvectionSolver2D *advecSolver;
};

#endif
