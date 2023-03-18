#ifndef __INS_LS_3D_H
#define __INS_LS_3D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "solvers/3d/advection_solver.h"

#include "dg_mesh/dg_mesh_3d.h"

class LevelSetSolver3D {
public:
  LevelSetSolver3D(DGMesh3D *m);
  LevelSetSolver3D(DGMesh3D *m, const std::string &filename);
  ~LevelSetSolver3D();

  void init();

  void setBCTypes(op_dat bc);
  void step(op_dat u, op_dat v, op_dat w, DG_FP dt);
  void getRhoMu(op_dat rho, op_dat mu);
  void getNormalsCurvature(op_dat nx, op_dat ny, op_dat nz, op_dat curv);
  void getDiracDelta(op_dat delta_x, op_dat delta_y, op_dat delta_z);

  DGMesh3D *mesh;

  op_dat s, s_modal, sampleX, sampleY, sampleZ, bc_types, smoothing;

  DG_FP alpha, order_width;
private:
  void sampleInterface();
  void reinitLS();
  // bool reinitNeeded();

  DG_FP h, epsilon, reinit_dt, reinit_width;
  int reinit_count;
  bool resuming;

  AdvectionSolver3D *advectionSolver;
};

#endif
