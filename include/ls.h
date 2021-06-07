#ifndef __INS_LS_H
#define __INS_LS_H

#include "op_seq.h"

#include "ins_data.h"

class LS {
public:
  LS(INSData *d);
  ~LS();

  void init();

  // void setVelField(op_dat u1, op_dat v1);
  // void step(double dt);

  INSData *data;
  op_dat u, v;
  op_dat s; //, s_bc;
  op_dat rk[3], rkQ;
  /*
  op_dat dsdx, dsdy, sign, dsldx, dsrdx, dsldy, dsrdy, dpldx, dprdx, dpldy, dprdy;
  op_dat sigmax, sigmay, sigmaFx, sigmaFy, diff, diffF, tau;
  */
private:
  /*
  void calc_sigma(double epsilon);
  void calc_diff(double epsilon);
  double get_residual();
  void reinit_ls();


  double h;
  double tol;
  int maxIter;
  double prev_res;
  */

  double *s_data;
  double *s_bc_data;
  double *rk_data[3];
  double *rkQ_data;
  /*
  double *dsdx_data;
  double *dsdy_data;
  double *sign_data;
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
  double *diff_data;
  double *diffF_data;
  double *tau_data;
  */
};

#endif
