#ifndef __DG_MATRIX_3_VEC_H
#define __DG_MATRIX_3_VEC_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

class Matrix3Vec {
public:
  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void set_bc_types(op_dat u_bc_ty, op_dat v_bc_ty, op_dat w_bc_ty) = 0;
  virtual void apply_bc(op_dat u_rhs, op_dat v_rhs, op_dat w_rhs, op_dat u_bc, op_dat v_bc, op_dat w_bc) = 0;
  virtual void mult(op_dat u_in, op_dat v_in, op_dat w_in, op_dat u_out, op_dat v_out, op_dat w_out) = 0;
};

#endif