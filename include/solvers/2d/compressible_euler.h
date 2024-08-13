#ifndef __INS_COMPRESSIBLE_EULER_2D_H
#define __INS_COMPRESSIBLE_EULER_2D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "dg_mesh/dg_mesh_2d.h"

class CompressibleEuler2D {
public:
  CompressibleEuler2D(DGMesh2D *m);
  ~CompressibleEuler2D();

  void init();
  void step();
  DG_FP get_time();
  void set_dt(const DG_FP _dt);
  
  op_dat Q[4];

private:
  void rhs(op_dat *inQ, op_dat *outQ);

  DGMesh2D *mesh;
  DG_FP h, dt, time;
  op_dat rk_wQ[4], rk_RHSQ[3][4];
};

#endif
