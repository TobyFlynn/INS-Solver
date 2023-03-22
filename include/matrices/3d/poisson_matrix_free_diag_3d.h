#ifndef __INS_POISSON_MATRIX_FREE_DIAG_3D_H
#define __INS_POISSON_MATRIX_FREE_DIAG_3D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_3d.h"
#include "../poisson_matrix_free_diag.h"

class PoissonMatrixFreeDiag3D : public PoissonMatrixFreeDiag {
public:
  PoissonMatrixFreeDiag3D(DGMesh3D *m);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;
  virtual void calc_mat_partial() override;

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  DGMesh3D *mesh;
  op_dat in_grad[3], tmp_npf[4], l[3];
};

#endif
