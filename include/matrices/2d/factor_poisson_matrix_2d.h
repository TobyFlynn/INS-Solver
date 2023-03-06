#ifndef __INS_FACTOR_POISSON_MATRIX_2D_H
#define __INS_FACTOR_POISSON_MATRIX_2D_H

#include "poisson_matrix_2d.h"

class FactorPoissonMatrix2D : public PoissonMatrix2D {
public:
  FactorPoissonMatrix2D(DGMesh2D *m);

  void set_factor(op_dat f);
  void set_surface(op_dat surf, const DG_FP al);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  op_dat factor, gFactor, cFactor, s, gS;
  DG_FP alpha;
};

#endif
