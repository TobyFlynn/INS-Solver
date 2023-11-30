#ifndef __INS_DIFFUSION_SOLVER_3D_H
#define __INS_DIFFUSION_SOLVER_3D_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"
#include "op_seq.h"

class DiffusionSolver3D {
public:
  DiffusionSolver3D(DGMesh3D *m);

  void step(op_dat val, op_dat vis);
  void set_dt();
  void set_dt(const DG_FP t);
  void set_dt(op_dat vis);
  DG_FP get_dt();
private:
  void rhs(op_dat val, op_dat vis, op_dat val_out);
  void rhs_over_int(op_dat val, op_dat vis, op_dat val_out);

  DGMesh3D *mesh;
  DG_FP dt;
  // bool over_int_diff;
};

#endif
