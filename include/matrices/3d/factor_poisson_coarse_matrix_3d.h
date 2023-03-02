#ifndef __INS_FACTOR_POISSON_COARSE_MATRIX_3D_H
#define __INS_FACTOR_POISSON_COARSE_MATRIX_3D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_3d.h"
#include "poisson_coarse_matrix_3d.h"

class FactorPoissonCoarseMatrix3D : public PoissonCoarseMatrix3D {
public:
  FactorPoissonCoarseMatrix3D(DGMesh3D *m);

  void set_factor(op_dat f);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  op_dat factor;
};

#endif
