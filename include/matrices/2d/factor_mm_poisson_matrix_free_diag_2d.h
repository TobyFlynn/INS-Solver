#ifndef __INS_MM_POISSON_MATRIX_FREE_DIAG_3D_H
#define __INS_MM_POISSON_MATRIX_FREE_DIAG_3D_H

#include "factor_poisson_matrix_free_diag_2d.h"

class FactorMMPoissonMatrixFreeDiag2D : public FactorPoissonMatrixFreeDiag2D {
public:
  FactorMMPoissonMatrixFreeDiag2D(DGMesh2D *m);

  void set_mm_factor(op_dat f);
  virtual void mult(op_dat in, op_dat out) override;
  virtual void calc_mat_partial() override;

private:
  void calc_mm();
  op_dat mm_factor;
};

#endif
