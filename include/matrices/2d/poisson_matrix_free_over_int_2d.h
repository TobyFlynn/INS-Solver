#ifndef __INS_POISSON_MATRIX_FREE_OVER_INT_2D_H
#define __INS_POISSON_MATRIX_FREE_OVER_INT_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "../poisson_matrix_free.h"

class PoissonMatrixFreeOverInt2D : public PoissonMatrixFree {
public:
  PoissonMatrixFreeOverInt2D(DGMesh2D *m);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;

protected:
  DGMesh2D *mesh;
  op_dat in_grad[2], gIn, gIn_grad[2], g_tmp[3], l[2];
};

#endif
