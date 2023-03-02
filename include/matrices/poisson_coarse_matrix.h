#ifndef __INS_POISSON_COARSE_MATRIX_H
#define __INS_POISSON_COARSE_MATRIX_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "petscvec.h"
#include "petscksp.h"

#include "dg_mesh/dg_mesh.h"
#include "poisson_matrix.h"

class PoissonCoarseMatrix : public PoissonMatrix {
public:
  virtual void mult(op_dat in, op_dat out) override;
  virtual void multJacobi(op_dat in, op_dat out) override;
  virtual int getUnknowns() override;

protected:
  virtual void set_glb_ind() override;
  virtual void setPETScMatrix() override;
};

#endif
