#ifndef __INS_FACTOR_POISSON_MATRIX_OVER_INT_2D_H
#define __INS_FACTOR_POISSON_MATRIX_OVER_INT_2D_H

#include "poisson_matrix_over_int_2d.h"

class FactorPoissonMatrixOverInt2D : public PoissonMatrixOverInt2D {
public:
  FactorPoissonMatrixOverInt2D(DGMesh2D *m);

  void set_factor(op_dat f);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  op_dat factor, gFactor;
};

#endif
