#ifndef __INS_POISSON_SEMI_MATRIX_FREE_2D_H
#define __INS_POISSON_SEMI_MATRIX_FREE_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "../poisson_semi_matrix_free.h"

class PoissonSemiMatrixFree2D : public PoissonSemiMatrixFree {
public:
  PoissonSemiMatrixFree2D(DGMesh2D *m);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;
  virtual void calc_mat_partial() override;

  // OP2 Dats
  op_dat h;
  op_dat glb_indBC;
  op_dat orderBC;

protected:
  virtual void calc_op1() override;
  virtual void calc_op2() override;
  virtual void calc_opbc() override;
  virtual void calc_glb_ind() override;

  DGMesh2D *mesh;
  op_dat in_grad[2], gIn, gIn_grad[2], g_tmp[3], l[2];
};

#endif
