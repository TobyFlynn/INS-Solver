#ifndef __INS_POISSON_SEMI_MATRIX_FREE_H
#define __INS_POISSON_SEMI_MATRIX_FREE_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "petscvec.h"
#include "petscksp.h"

#include "dg_mesh/dg_mesh.h"
#include "poisson_matrix.h"

class PoissonSemiMatrixFree : public PoissonMatrix {
public:
  virtual void calc_mat_partial() = 0;
  virtual void calc_mat() override;
  virtual void mult(op_dat in, op_dat out) override;
  virtual void multJacobi(op_dat in, op_dat out) override;
  virtual bool getPETScMat(Mat** mat) override;

protected:
  virtual void setPETScMatrix() override;
};

#endif
