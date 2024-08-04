#ifndef __DG_VISCOUS_SOLVER_H
#define __DG_VISCOUS_SOLVER_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "matrix_2_vec.h"
#include "dg_mesh/dg_mesh_2d.h"

class ViscousSolver {
public:
  enum Preconditioners {
    NONE, INV_MASS, FACTOR_INV_MASS, RECP_FACTOR_DAT_INV_MASS, JACOBI
  };
  ViscousSolver(DGMesh2D *m);
  virtual void set_matrix(Matrix2Vec *mat);
  void set_bcs(op_dat u_bcs, op_dat v_bcs);
  void set_nullspace(bool ns);
  virtual bool solve(op_dat u_rhs, op_dat v_rhs, op_dat u_ans, op_dat v_ans);
  virtual void set_tol_and_iter(const double rtol, const double atol, const int maxiter);
  virtual void set_preconditioner(Preconditioners p);
  void set_inv_mass_factor(DG_FP f);
  void set_inv_mass_recp_factor(op_dat f);

protected:
  Matrix2Vec *matrix;
  bool nullspace, zero_input;
  op_dat u_bc, v_bc, inv_mass_recp_factor;
  DGMesh2D *mesh;
  Preconditioners preconditioner;
  int max_iter;
  DG_FP abs_tol, rel_tol;
  DG_FP inv_mass_factor;

  virtual bool conjugate_gradient(op_dat u_rhs, op_dat v_rhs, op_dat u_ans, op_dat v_ans);
  virtual bool preconditioned_conjugate_gradient(op_dat u_rhs, op_dat v_rhs, op_dat u_ans, op_dat v_ans);

  // Helper functions
  void update_ans(DG_FP alpha, op_dat p0, op_dat p1, op_dat u, op_dat v);
  DG_FP update_res(DG_FP alpha, op_dat Ap0, op_dat Ap1, op_dat r0, op_dat r1);
  DG_FP calc_res_explicit(op_dat u_rhs, op_dat v_rhs, op_dat u_ans, op_dat v_ans, op_dat r0, op_dat r1);
  void update_p(DG_FP beta, op_dat r0, op_dat r1, op_dat p0, op_dat p1);

  // Preconditioner functions
  virtual void precondition(op_dat u_res, op_dat v_res, op_dat u, op_dat v);
  virtual void pre_inv_mass(op_dat u_res, op_dat v_res, op_dat u, op_dat v);
  virtual void pre_factor_inv_mass(op_dat u_res, op_dat v_res, op_dat u, op_dat v);
  virtual void pre_recp_dat_factor_inv_mass(op_dat u_res, op_dat v_res, op_dat u, op_dat v);
  virtual void pre_jacobi(op_dat u_res, op_dat v_res, op_dat u, op_dat v);
};

#endif