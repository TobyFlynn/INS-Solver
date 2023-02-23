#ifndef __INS_FACTOR_POISSON_SEMI_MATRIX_FREE_3D_H
#define __INS_FACTOR_POISSON_SEMI_MATRIX_FREE_3D_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_3d.h"
#include "factor_poisson_matrix_3d.h"

class FactorPoissonSemiMatrixFree3D : public FactorPoissonMatrix3D {
public:
  FactorPoissonSemiMatrixFree3D(DGMesh3D *m, bool init_mat_dats = true);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void mult(op_dat in, op_dat out) override;

protected:
  op_dat in_grad[3], tmp_npf[4], l[3];
};

#endif
