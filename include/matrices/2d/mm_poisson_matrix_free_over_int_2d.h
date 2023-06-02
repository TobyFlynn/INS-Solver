#ifndef __INS_MM_POISSON_MATRIX_FREE_OVER_INT_2D_H
#define __INS_MM_POISSON_MATRIX_FREE_OVER_INT_2D_H

#include "poisson_matrix_free_over_int_2d.h"

class MMPoissonMatrixFreeOverInt2D : public PoissonMatrixFreeOverInt2D {
public:
  MMPoissonMatrixFreeOverInt2D(DGMesh2D *m);

  void set_factor(DG_FP f);
  DG_FP get_factor();
  virtual void mult(op_dat in, op_dat out) override;

private:
  DG_FP factor;
};

#endif
