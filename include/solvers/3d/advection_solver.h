#ifndef __INS_ADVECTION_SOLVER_3D_H
#define __INS_ADVECTION_SOLVER_3D_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"
#include "op_seq.h"

class AdvectionSolver3D {
public:
  AdvectionSolver3D(DGMesh3D *m);

  void step(op_dat val, op_dat u, op_dat v, op_dat w);
  void set_dt();
  void set_dt(const DG_FP t);
  void set_bc_types(op_dat bc);

protected:
  virtual void bc_kernel(op_dat val, op_dat u, op_dat v, op_dat w,
                         op_dat out) = 0;
  virtual void bc_kernel_oi(op_dat val, op_dat u, op_dat v, op_dat w,
                            op_dat flux) = 0;

  DGMesh3D *mesh;

private:
  void rhs(op_dat val, op_dat u, op_dat v, op_dat w, op_dat val_out);
  void rhs_over_int(op_dat val, op_dat u, op_dat v, op_dat w, op_dat val_out);

  bool over_int_advec;
  DG_FP dt;
  op_dat bc_types;
};

#endif
