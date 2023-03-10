#ifndef __INS_CUB_FACTOR_MM_POISSON_MATRIX_2D_H
#define __INS_CUB_FACTOR_MM_POISSON_MATRIX_2D_H

#include "cub_factor_poisson_matrix_2d.h"

class CubFactorMMPoissonMatrix2D : public CubFactorPoissonMatrix2D {
public:
  CubFactorMMPoissonMatrix2D(DGMesh2D *m);

  virtual void calc_mat() override;
  void set_mm_factor(op_dat f);
private:
  void calc_mm();
  op_dat mm_factor;
};

#endif
