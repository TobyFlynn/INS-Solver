#ifndef __INS_DIFFUSION_SOLVER_2D_H
#define __INS_DIFFUSION_SOLVER_2D_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_2d.h"
#include "op_seq.h"

class DiffusionSolver2D {
public:
  DiffusionSolver2D(DGMesh2D *m);

  void step(op_dat val, op_dat vis);
  void set_dt();
  void set_dt(const DG_FP t);
private:
  void rhs(op_dat val, op_dat vis, op_dat val_out);
  void rhs_over_int(op_dat val, op_dat vis, op_dat val_out);

  DGMesh2D *mesh;
  DG_FP dt;
  bool over_int_diff;
};

#endif
