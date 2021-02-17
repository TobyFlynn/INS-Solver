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
  Poisson(INSData *data);
  ~Poisson();
private:
  // OP2 Dats
  op_dat pTau, pExRHS[2];
  // Pointers to private memory
  double *pTau_data;
  double *pExRHS_data[2];
};

#endif
