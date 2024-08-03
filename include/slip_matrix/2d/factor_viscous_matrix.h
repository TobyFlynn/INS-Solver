#ifndef __DG_FACTOR_VISCOUS_MATRIX_2D_H
#define __DG_FACTOR_VISCOUS_MATRIX_2D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "dg_mesh/dg_mesh_2d.h"
#include "matrix_2_vec.h"

class FactorViscousMatrix2D : public Matrix2Vec {
public:
  FactorViscousMatrix2D(DGMesh2D *m);

  // op_dat bc_types - 0 for Dirichlet, 1 for Neumann
  virtual void set_bc_types(op_dat u_bc_ty, op_dat v_bc_ty) override;
  virtual void apply_bc(op_dat u_rhs, op_dat v_rhs, op_dat u_bc, op_dat v_bc) override;
  virtual void mult(op_dat u_in, op_dat v_in, op_dat u_out, op_dat v_out) override;
  void set_factor(op_dat f);
  void set_mm_factor(op_dat f);
  // virtual void mult_sp(op_dat in, op_dat out);
  // virtual void multJacobi(op_dat in, op_dat out);
  // virtual void multJacobi_sp(op_dat in, op_dat out);
  // virtual bool getPETScMat(Mat** mat);
  // virtual DG_MAT_IND_TYPE getUnknowns();

  // op_dat op1, op2[2], opbc, glb_ind, glb_indL, glb_indR;
  // int unknowns;
protected:
  // virtual void set_glb_ind();
  // virtual void calc_glb_ind() = 0;
  // virtual void setPETScMatrix();
  virtual void mat_free_pre_compute_tau();

  DGMesh2D *mesh;
  op_dat u_bc_types, v_bc_types, mat_free_tau_c, factor, mm_factor;
};

#endif