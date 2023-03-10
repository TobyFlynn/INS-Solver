#ifndef __INS_CUB_FACTOR_POISSON_MATRIX_2D_H
#define __INS_CUB_FACTOR_POISSON_MATRIX_2D_H

#include "cub_poisson_matrix_2d.h"

class CubFactorPoissonMatrix2D : public CubPoissonMatrix2D {
public:
  CubFactorPoissonMatrix2D(DGMesh2D *m);

  void set_factor(op_dat f);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  op_dat factor, gFactor, cFactor;
};

#endif
