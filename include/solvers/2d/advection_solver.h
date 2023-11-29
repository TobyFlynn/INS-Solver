#ifndef __INS_ADVECTION_SOLVER_2D_H
#define __INS_ADVECTION_SOLVER_2D_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_2d.h"
#include "op_seq.h"

class AdvectionSolver2D {
public:
  AdvectionSolver2D(DGMesh2D *m);

  void step(op_dat val, op_dat u, op_dat v);
  void set_dt();
  void set_dt(const DG_FP t);

protected:
  virtual void bc_kernel(op_dat val, op_dat u, op_dat v, op_dat out) = 0;
  virtual void bc_kernel_oi(op_dat val, op_dat u, op_dat v, op_dat flux) = 0;

  DGMesh2D *mesh;

private:
  void rhs(op_dat val, op_dat u, op_dat v, op_dat val_out);
  void rhs_over_int(op_dat val, op_dat u, op_dat v, op_dat val_out);

  DG_FP dt;
  bool over_int_advec;
};

#endif
