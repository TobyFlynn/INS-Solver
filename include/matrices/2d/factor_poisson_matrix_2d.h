#ifndef __INS_FACTOR_POISSON_MATRIX_2D_H
#define __INS_FACTOR_POISSON_MATRIX_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "poisson_matrix_2d.h"

class FactorPoissonMatrix2D : public PoissonMatrix2D {
public:
  FactorPoissonMatrix2D(DGMesh2D *m);

  void set_factor(op_dat f);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc(op_dat bc_types) override;

  op_dat factor, gFactor, cFactor;
};

#endif
