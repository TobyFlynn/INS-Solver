#ifndef __INS_FACTOR_COARSE_POISSON_MATRIX_OVER_INT_2D_H
#define __INS_FACTOR_COARSE_POISSON_MATRIX_OVER_INT_2D_H

#include "poisson_coarse_matrix_over_int_2d.h"

class FactorPoissonCoarseMatrixOverInt2D : public PoissonCoarseMatrixOverInt2D {
public:
  FactorPoissonCoarseMatrixOverInt2D(DGMesh2D *m);

  void set_factor(op_dat f);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  op_dat factor, gFactor;
};

#endif
