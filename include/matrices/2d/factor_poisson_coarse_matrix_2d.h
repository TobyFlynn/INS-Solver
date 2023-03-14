#ifndef __INS_FACTOR_COARSE_POISSON_MATRIX_2D_H
#define __INS_FACTOR_COARSE_POISSON_MATRIX_2D_H

#include "poisson_coarse_matrix_2d.h"

class FactorPoissonCoarseMatrix2D : public PoissonCoarseMatrix2D {
public:
  FactorPoissonCoarseMatrix2D(DGMesh2D *m);

  void set_factor(op_dat f);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  op_dat factor, gFactor;
};

#endif
