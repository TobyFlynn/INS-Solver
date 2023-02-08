#ifndef __INS_FACTOR_MM_POISSON_MATRIX_2D_H
#define __INS_FACTOR_MM_POISSON_MATRIX_2D_H

#include "factor_poisson_matrix_2d.h"

class FactorMMPoissonMatrix2D : public FactorPoissonMatrix2D {
public:
  FactorMMPoissonMatrix2D(DGMesh2D *m);

  virtual void calc_mat(op_dat bc_types) override;
  void set_mm_factor(op_dat f);
private:
  void calc_mm();
  op_dat mm_factor;
};

#endif
