#ifndef __INS_MM_POISSON_MATRIX_2D_H
#define __INS_MM_POISSON_MATRIX_2D_H

#include "poisson_matrix_2d.h"

class MMPoissonMatrix2D : public PoissonMatrix2D {
public:
  MMPoissonMatrix2D(DGMesh2D *m);

  virtual void calc_mat() override;
  void set_factor(DG_FP f);
  DG_FP get_factor();

private:
  void calc_mm();

  DG_FP factor;
};

#endif
