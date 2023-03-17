#ifndef __INS_FACTOR_POISSON_MATRIX_FREE_2D_H
#define __INS_FACTOR_POISSON_MATRIX_FREE_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "poisson_matrix_free_2d.h"

class FactorPoissonMatrixFree2D : public PoissonMatrixFree2D {
public:
  FactorPoissonMatrixFree2D(DGMesh2D *m);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;
  void set_factor(op_dat f);

protected:
  op_dat factor, gFactor;
};

#endif
