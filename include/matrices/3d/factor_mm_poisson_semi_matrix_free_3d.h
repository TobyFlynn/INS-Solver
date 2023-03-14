#ifndef __INS_MM_POISSON_MATRIX_FREE_3D_H
#define __INS_MM_POISSON_MATRIX_FREE_3D_H

#include "factor_poisson_semi_matrix_free_3d.h"

class FactorMMPoissonSemiMatrixFree3D : public FactorPoissonSemiMatrixFree3D {
public:
  FactorMMPoissonSemiMatrixFree3D(DGMesh3D *m);

  void set_mm_factor(op_dat f);
  virtual void mult(op_dat in, op_dat out) override;
  virtual void calc_mat_partial() override;

private:
  void calc_mm();
  op_dat mm_factor;
};

#endif
