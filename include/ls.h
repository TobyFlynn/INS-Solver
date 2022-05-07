#ifndef __INS_LS_H
#define __INS_LS_H

#include "op_seq.h"

#include "ins_data.h"

#include "dg_global_constants.h"

class LS {
public:
  LS(DGMesh *m, INSData *d);
  ~LS();

  void init();

  void setVelField(op_dat u1, op_dat v1);
  void step(double dt);

  DGMesh *mesh;
  INSData *data;

  op_dat u, v;
  op_dat s, step_s, nx, ny, curv, diracDelta;
  op_dat rk[3], rkQ;
  op_dat F, G, dFdr, dFds, dGdr, dGds, gInput, gU, gV, nFlux, exAdvec;

  op_dat dsdx, dsdy;

  double alpha, order_width;
private:
  void advec_step(op_dat input, op_dat output);
  void reinit_ls();
  bool reinit_needed();
  void update_values();

  double h, epsilon, reinit_dt;
  int numSteps;

  double *s_data, *step_s_data, *nx_data, *ny_data, *curv_data, *diracDelta_data;
  double *gInput_data;
};

#endif
