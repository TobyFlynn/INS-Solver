#ifndef __INS_FACTOR_MM_POISSON_MATRIX_OVER_INT_2D_H
#define __INS_FACTOR_MM_POISSON_MATRIX_OVER_INT_2D_H

#include "factor_poisson_matrix_over_int_2d.h"

class FactorMMPoissonMatrixOverInt2D : public FactorPoissonMatrixOverInt2D {
public:
  FactorMMPoissonMatrixOverInt2D(DGMesh2D *m);

  virtual void calc_mat() override;
  void set_mm_factor(op_dat f);
private:
  void calc_mm();
  op_dat mm_factor;
};

#endif
