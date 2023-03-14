#ifndef __INS_CUB_FACTOR_POISSON_MATRIX_2D_H
#define __INS_CUB_FACTOR_POISSON_MATRIX_2D_H

#include "factor_poisson_matrix_2d.h"

class CubFactorPoissonMatrix2D : public FactorPoissonMatrix2D {
public:
  CubFactorPoissonMatrix2D(DGMesh2D *m);

protected:
  virtual void calc_op1() override;

  op_dat cFactor;
};

#endif
