#ifndef __INS_POISSON_MATRIX_H
#define __INS_POISSON_MATRIX_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "petscvec.h"
#include "petscksp.h"

#include "dg_mesh/dg_mesh.h"

class PoissonMatrix {
public:
  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void calc_mat() = 0;
  virtual void set_bc_types(op_dat bc_ty);
  virtual void apply_bc(op_dat rhs, op_dat bc) = 0;
  virtual void mult(op_dat in, op_dat out);
  virtual void multJacobi(op_dat in, op_dat out);
  virtual bool getPETScMat(Mat** mat);

  op_dat op1, op2[2], opbc, glb_ind, glb_indL, glb_indR, orderL, orderR;
  int unknowns;
protected:
  virtual void calc_op1() = 0;
  virtual void calc_op2() = 0;
  virtual void calc_opbc() = 0;
  virtual void set_glb_ind();
  virtual void calc_glb_ind() = 0;
  virtual void setPETScMatrix();

  DGMesh *_mesh;

  Mat pMat;
  bool petscMatInit;
  bool petscMatResetRequired;

  op_dat bc_types;
};

#endif
