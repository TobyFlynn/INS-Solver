#ifndef __INS_MM_POISSON_MATRIX_3D_H
#define __INS_MM_POISSON_MATRIX_3D_H

#include "poisson_matrix_3d.h"

class MMPoissonMatrix3D : public PoissonMatrix3D {
public:
  MMPoissonMatrix3D(DGMesh3D *m);

  virtual void calc_mat() override;
  void set_factor(double f);
  double get_factor();

private:
  void calc_mm();

  double factor;
};

#endif
