#ifndef __INS_MM_POISSON_MATRIX_FREE_3D_H
#define __INS_MM_POISSON_MATRIX_FREE_3D_H

#include "factor_poisson_matrix_free_3d.h"

class FactorMMPoissonMatrixFree3D : public FactorPoissonMatrixFree3D {
public:
  FactorMMPoissonMatrixFree3D(DGMesh3D *m);

  void set_mm_factor(op_dat f);
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;

private:
  op_dat mm_factor;
};

#endif
