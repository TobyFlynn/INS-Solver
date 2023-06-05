#ifndef __INS_POISSON_MATRIX_FREE_2D_H
#define __INS_POISSON_MATRIX_FREE_2D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "../poisson_matrix_free.h"
#include "poisson_matrix_free_mult_2d.h"

class PoissonMatrixFree2D : public PoissonMatrixFree, public PoissonMatrixFreeMult2D {
public:
  PoissonMatrixFree2D(DGMesh2D *m);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void set_bc_types(op_dat bc_ty) override;
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;
};

#endif
