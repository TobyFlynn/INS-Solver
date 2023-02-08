#ifndef __INS_POISSON_MATRIX_H
#define __INS_POISSON_MATRIX_H

#include "op_seq.h"

#include "petscvec.h"
#include "petscksp.h"

class PoissonMatrix {
public:
  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void calc_mat(op_dat bc_types) = 0;
  virtual void apply_bc(op_dat rhs, op_dat bc) = 0;
  virtual void mult(op_dat in, op_dat out) = 0;
  virtual void multJacobi(op_dat in, op_dat out) = 0;
  virtual bool getPETScMat(Mat** mat) = 0;

  op_dat op1, op2[2], opbc, glb_ind, glb_indL, glb_indR;
protected:
  virtual void calc_op1() = 0;
  virtual void calc_op2() = 0;
  virtual void calc_opbc(op_dat bc_types) = 0;
  virtual void set_glb_ind() = 0;
  virtual void calc_glb_ind() = 0;

  Mat pMat;
  bool petscMatInit;
  bool petscMatResetRequired;
};

#endif
