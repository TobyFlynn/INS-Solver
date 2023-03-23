#ifndef __INS_POISSON_SEMI_MATRIX_FREE_3D_H
#define __INS_POISSON_SEMI_MATRIX_FREE_3D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_3d.h"
#include "../poisson_semi_matrix_free.h"

class PoissonSemiMatrixFree3D : public PoissonSemiMatrixFree {
public:
  PoissonSemiMatrixFree3D(DGMesh3D *m, bool alloc_tmp_dats = true);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;
  virtual void calc_mat_partial() override;
  virtual void set_tmp_dats(op_dat np0, op_dat np1, op_dat np2, op_dat npf0,
                            op_dat npf1, op_dat npf2, op_dat npf3);

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;

  DGMesh3D *mesh;
  op_dat in_grad[3], tmp_npf[4], l[3];
};

#endif
