#ifndef __INS_ADVECTION_SOLVER_H
#define __INS_ADVECTION_SOLVER_H

#include "dg_mesh/dg_mesh_2d.h"
#include "op_seq.h"

class AdvectionSolver2D {
public:
  AdvectionSolver2D(DGMesh2D *m);

  void step(op_dat val, op_dat u, op_dat v);
  void set_dt();
  void set_dt(const double t);
private:
  void rhs(op_dat val, op_dat u, op_dat v, op_dat val_out);

  DGMesh2D *mesh;
  double dt;
  op_dat f, g, flux, rk[3], rkQ, gVal, gU, gV;
};

#endif
