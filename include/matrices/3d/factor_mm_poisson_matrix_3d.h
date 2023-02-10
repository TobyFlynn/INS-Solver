#ifndef __INS_FACTOR_MM_POISSON_MATRIX_3D_H
#define __INS_FACTOR_MM_POISSON_MATRIX_3D_H

#include "factor_poisson_matrix_3d.h"

class FactorMMPoissonMatrix3D : public FactorPoissonMatrix3D {
public:
  FactorMMPoissonMatrix3D(DGMesh3D *m);

  virtual void calc_mat() override;
  void set_mm_factor(op_dat f);
private:
  void calc_mm();
  op_dat mm_factor;
};

#endif
