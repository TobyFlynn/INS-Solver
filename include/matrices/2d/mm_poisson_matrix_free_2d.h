#ifndef __INS_MM_POISSON_MATRIX_FREE_2D_H
#define __INS_MM_POISSON_MATRIX_FREE_2D_H

#include "poisson_matrix_free_2d.h"

class MMPoissonMatrixFree2D : public PoissonMatrixFree2D {
public:
  MMPoissonMatrixFree2D(DGMesh2D *m);

  void set_factor(DG_FP f);
  DG_FP get_factor();
  virtual void mult(op_dat in, op_dat out) override;

private:
  DG_FP factor;
};

#endif
