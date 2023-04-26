#ifndef __INS_FACTOR_POISSON_SEMI_MATRIX_FREE_3D_H
#define __INS_FACTOR_POISSON_SEMI_MATRIX_FREE_3D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_3d.h"
#include "../poisson_semi_matrix_free.h"
#include "factor_poisson_matrix_free_mult_3d.h"

class FactorPoissonSemiMatrixFree3D : public PoissonSemiMatrixFree, public FactorPoissonMatrixFreeMult3D {
public:
  FactorPoissonSemiMatrixFree3D(DGMesh3D *m);

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
