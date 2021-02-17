#ifndef __POISSON_H
#define __POISSON_H

#include "op_seq.h"
#include "ins_data.h"

extern double ones[15];
extern double r[15];
extern double s[15];
extern double Dr[15 * 15];
extern double Ds[15 * 15];
extern double Drw[15 * 15];
extern double Dsw[15 * 15];
extern double LIFT[15 * 15];
extern double MASS[15 * 15];
extern int FMASK[15];

class Poisson {
public:
  Poisson(INSData *nsData);
  ~Poisson();
  // OP2 Dats
  op_dat pTau, pExRHS[2], pU, pDu, pFluxXu, pFluxYu, pDuDx, pDuDy, pFluxQ, pDivQ, pRHS;
private:
  void rhs(const double *u, double *rhs);
  INSData *data;
  // Pointers to private memory
  double *pTau_data;
  double *pExRHS_data[2];
  double *pU_data;
  double *pDu_data;
  double *pFluxXu_data;
  double *pFluxYu_data;
  double *pDuDx_data;
  double *pDuDy_data;
  double *pFluxQ_data;
  double *pDivQ_data;
  double *pRHS_data;
};

#endif
