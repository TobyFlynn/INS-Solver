#ifndef __INS_ADVEC_DIFF_SOLVER_2D_H
#define __INS_ADVEC_DIFF_SOLVER_2D_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_2d.h"
#include "op_seq.h"

class AdvecDiffSolver2D {
public:
  AdvecDiffSolver2D(DGMesh2D *m);

  void step(op_dat val, op_dat u, op_dat v, op_dat vis, DG_FP time);
  void step(op_dat val, op_dat u, op_dat v, op_dat vis);
  void set_dt();
  void set_dt(const DG_FP t);
  void set_dt(op_dat u, op_dat v, op_dat vis);
  DG_FP get_dt();
private:
  void rhs(op_dat val, op_dat u, op_dat v, op_dat vis, op_dat val_out);
  void advec(op_dat val, op_dat u, op_dat v, op_dat val_out);
  void diff(op_dat val, op_dat vis, op_dat val_out);
  void rhs_over_int(op_dat val, op_dat vis, op_dat val_out);

  DGMesh2D *mesh;
  DG_FP dt;
  bool over_int_diff;
};

#endif
