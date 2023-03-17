#ifndef __INS_FACTOR_POISSON_SEMI_MATRIX_FREE_2D_H
#define __INS_FACTOR_POISSON_SEMI_MATRIX_FREE_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "poisson_semi_matrix_free_2d.h"

class FactorPoissonSemiMatrixFree2D : public PoissonSemiMatrixFree2D {
public:
  FactorPoissonSemiMatrixFree2D(DGMesh2D *m);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;
  void set_factor(op_dat f);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  op_dat factor, gFactor;
};

#endif
