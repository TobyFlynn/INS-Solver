#ifndef __INS_POISSON_MATRIX_FREE_MULT_2D_H
#define __INS_POISSON_MATRIX_FREE_MULT_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"

class PoissonMatrixFreeMult2D {
public:
  PoissonMatrixFreeMult2D(DGMesh2D *m);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void mat_free_set_bc_types(op_dat bc_ty);
  virtual void mat_free_apply_bc(op_dat rhs, op_dat bc);
  virtual void mat_free_mult(op_dat in, op_dat out);

protected:
  DGMesh2D *mesh;
  op_dat mat_free_bcs;
};

#endif
