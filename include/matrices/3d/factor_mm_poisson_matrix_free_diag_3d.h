#ifndef __INS_MM_POISSON_MATRIX_FREE_DIAG_3D_H
#define __INS_MM_POISSON_MATRIX_FREE_DIAG_3D_H

#include "factor_poisson_matrix_free_diag_3d.h"

class FactorMMPoissonMatrixFreeDiag3D : public FactorPoissonMatrixFreeDiag3D {
public:
  FactorMMPoissonMatrixFreeDiag3D(DGMesh3D *m, bool alloc_tmp_dats = true);

  void set_mm_factor(op_dat f);
  virtual void mult(op_dat in, op_dat out) override;
  virtual void calc_mat_partial() override;

private:
  void calc_mm();
  op_dat mm_factor;
};

#endif
