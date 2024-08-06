#ifndef __DG_FACTOR_VISCOUS_MATRIX_3D_H
#define __DG_FACTOR_VISCOUS_MATRIX_3D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "dg_mesh/dg_mesh_3d.h"
#include "matrix_3_vec.h"

class FactorViscousMatrix3D : public Matrix3Vec {
public:
  FactorViscousMatrix3D(DGMesh3D *m, bool calc_diagonal = false, bool calc_inv_block_diagonal = false);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void set_bc_types(op_dat u_bc_ty, op_dat v_bc_ty, op_dat w_bc_ty) override;
  virtual void apply_bc(op_dat u_rhs, op_dat v_rhs, op_dat w_rhs, op_dat u_bc, op_dat v_bc, op_dat w_bc) override;
  virtual void mult(op_dat u_in, op_dat v_in, op_dat w_in, op_dat u_out, op_dat v_out, op_dat w_out) override;
  void set_factor(op_dat f);
  void set_mm_factor(op_dat f);
  void calc_diag();
  void calc_inv_block_diag();

  op_dat u_diag, v_diag, w_diag, u_inv_block_diag, v_inv_block_diag, w_inv_block_diag;
protected:
  virtual void mat_free_pre_compute_tau();

  DGMesh3D *mesh;
  op_dat u_bc_types, v_bc_types, w_bc_types, mat_free_tau_c, factor, mm_factor;
};

#endif