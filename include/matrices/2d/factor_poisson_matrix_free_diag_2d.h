#ifndef __INS_FACTOR_POISSON_MATRIX_FREE_DIAG_2D_H
#define __INS_FACTOR_POISSON_MATRIX_FREE_DIAG_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "../poisson_matrix_free_diag.h"
#include "factor_poisson_matrix_free_mult_2d.h"

class FactorPoissonMatrixFreeDiag2D : public PoissonMatrixFreeDiag, public FactorPoissonMatrixFreeMult2D {
public:
  FactorPoissonMatrixFreeDiag2D(DGMesh2D *m);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void set_bc_types(op_dat bc_ty) override;
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;
  virtual void calc_mat_partial() override;
  void set_factor(op_dat f);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  op_dat factor;
};

#endif
