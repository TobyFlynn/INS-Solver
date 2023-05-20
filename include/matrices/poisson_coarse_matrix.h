#ifndef __INS_POISSON_COARSE_MATRIX_H
#define __INS_POISSON_COARSE_MATRIX_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "petscvec.h"
#include "petscksp.h"

#include "dg_mesh/dg_mesh.h"
#include "poisson_matrix.h"

#ifdef INS_CUDA
#include <amgx_c.h>
#endif

class PoissonCoarseMatrix : public PoissonMatrix {
public:
  virtual void mult(op_dat in, op_dat out) override;
  virtual void multJacobi(op_dat in, op_dat out) override;
  virtual int getUnknowns() override;

  #ifdef INS_CUDA
  bool getAmgXMat(AMGX_matrix_handle** mat);
  #endif

protected:
  virtual void set_glb_ind() override;
  virtual void setPETScMatrix() override;
  #ifdef INS_CUDA
  virtual void setAmgXMatrix();

  AMGX_matrix_handle amgx_mat;
  bool amgx_mat_init = false;
  #endif
};

#endif
