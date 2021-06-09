#ifndef __INS_LS_H
#define __INS_LS_H

#include "op_seq.h"

#include "ins_data.h"

class LS {
public:
  LS(INSData *d, CubatureData *c, GaussData *g);
  ~LS();

  void init();

  void setVelField(op_dat u1, op_dat v1);
  void step(double dt);

  INSData *data;
  CubatureData *cData;
  GaussData *gData;

  op_dat u, v;
  op_dat s; //, s_bc;
  op_dat rk[3], rkQ;
  op_dat F, G, dFdr, dFds, dGdr, dGds, nFlux, exAdvec;

  op_dat dsdx, dsdy, sign, gS, dsldx, dsrdx, dsldy, dsrdy, dpldx, dprdx, dpldy, dprdy;
  op_dat sigmax, sigmay, sigmaFx, sigmaFy, gSigmax, gSigmay, diff, diffF;
private:
  void advec_step(op_dat input, op_dat output);

  // void calc_sigma(double epsilon);
  void calc_diff(double epsilon);
  // double get_residual();
  void reinit_ls();

  double h;
  double alpha;
  double epsilon;
  double reinit_dt;
  int numSteps;
  int counter;

  double *s_data;
  double *s_bc_data;
  double *rk_data[3];
  double *rkQ_data;

  double *F_data;
  double *G_data;
  double *dFdr_data;
  double *dFds_data;
  double *dGdr_data;
  double *dGds_data;
  double *nFlux_data;
  double *exAdvec_data;

  double *dsdx_data;
  double *dsdy_data;
  double *sign_data;
  double *gS_data;
  double *dsldx_data;
  double *dsrdx_data;
  double *dsldy_data;
  double *dsrdy_data;
  double *dpldx_data;
  double *dprdx_data;
  double *dpldy_data;
  double *dprdy_data;

  double *sigmax_data;
  double *sigmay_data;
  double *sigmaFx_data;
  double *sigmaFy_data;
  double *gSigmax_data;
  double *gSigmay_data;
  double *diff_data;
  double *diffF_data;
};

#endif
