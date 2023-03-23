#ifndef __INS_MM_POISSON_MATRIX_FREE_3D_H
#define __INS_MM_POISSON_MATRIX_FREE_3D_H

#include "poisson_matrix_free_3d.h"

class MMPoissonMatrixFree3D : public PoissonMatrixFree3D {
public:
  MMPoissonMatrixFree3D(DGMesh3D *m, bool alloc_tmp_dats = true);

  void set_factor(DG_FP f);
  DG_FP get_factor();
  virtual void mult(op_dat in, op_dat out) override;

private:
  DG_FP factor;
};

#endif
