#ifndef __INS_MM_POISSON_MATRIX_FREE_3D_H
#define __INS_MM_POISSON_MATRIX_FREE_3D_H

#include "poisson_matrix_3d.h"

class MMPoissonMatrixFree3D : public PoissonMatrix3D {
public:
  MMPoissonMatrixFree3D(DGMesh3D *m);

  virtual void calc_mat() override;
  void set_factor(DG_FP f);
  DG_FP get_factor();
  virtual void mult(op_dat in, op_dat out) override;

private:
  void calc_mm();

  DG_FP factor;
};

#endif
