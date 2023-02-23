#ifndef __INS_POISSON_SEMI_MATRIX_FREE_3D_H
#define __INS_POISSON_SEMI_MATRIX_FREE_3D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_3d.h"
#include "poisson_matrix_3d.h"

class PoissonSemiMatrixFree3D : public PoissonMatrix3D {
public:
  PoissonSemiMatrixFree3D(DGMesh3D *m, bool init_mat_dats = true);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void apply_bc(op_dat rhs, op_dat bc) override;
  virtual void mult(op_dat in, op_dat out) override;

protected:
  op_dat in_grad[3], tmp_npf[4], l[3];
};

#endif
