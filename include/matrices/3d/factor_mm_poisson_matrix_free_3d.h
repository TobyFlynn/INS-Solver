#ifndef __INS_MM_POISSON_MATRIX_FREE_3D_H
#define __INS_MM_POISSON_MATRIX_FREE_3D_H

#include "factor_poisson_matrix_free_3d.h"

class FactorMMPoissonMatrixFree3D : public FactorPoissonMatrixFree3D {
public:
  FactorMMPoissonMatrixFree3D(DGMesh3D *m, bool alloc_tmp_dats = true);

  void set_mm_factor(op_dat f);
  virtual void mult(op_dat in, op_dat out) override;

private:
  op_dat mm_factor;
};

#endif
