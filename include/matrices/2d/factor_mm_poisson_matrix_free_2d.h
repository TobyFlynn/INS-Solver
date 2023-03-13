#ifndef __INS_MM_POISSON_MATRIX_FREE_2D_H
#define __INS_MM_POISSON_MATRIX_FREE_2D_H

#include "factor_poisson_matrix_free_2d.h"

class FactorMMPoissonMatrixFree2D : public FactorPoissonMatrixFree2D {
public:
  FactorMMPoissonMatrixFree2D(DGMesh2D *m);

  void set_mm_factor(op_dat f);
  virtual void mult(op_dat in, op_dat out) override;

private:
  op_dat mm_factor;
};

#endif
