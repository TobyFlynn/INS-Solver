#ifndef __INS_FACTOR_POISSON_MATRIX_3D_H
#define __INS_FACTOR_POISSON_MATRIX_3D_H

#include "poisson_matrix_3d.h"

class FactorPoissonMatrix3D : public PoissonMatrix3D {
public:
  FactorPoissonMatrix3D(DGMesh3D *m, bool init_mat_dats = true);

  void set_factor(op_dat f);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  op_dat factor;
};

#endif
