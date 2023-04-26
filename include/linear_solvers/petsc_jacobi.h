#ifndef __PETSC_JACOBI_H
#define __PETSC_JACOBI_H

#include "op_seq.h"
#include "linear_solver.h"
#include "petscvec.h"
#include "petscksp.h"
#include "dg_mesh/dg_mesh.h"
#include "matrices/poisson_matrix_free_diag.h"

class PETScJacobiSolver : public LinearSolver {
public:
  PETScJacobiSolver(DGMesh *m);
  ~PETScJacobiSolver();

  bool solve(op_dat rhs, op_dat ans) override;

  void calc_rhs(const DG_FP *in_d, DG_FP *out_d);
  void precond(const DG_FP *in_d, DG_FP *out_d);

private:
  void create_shell_mat();
  void set_shell_pc(PC pc);

  DGMesh *mesh;
  PoissonMatrixFreeDiag *diagMat;
  KSP ksp;
  bool pMatInit, dat_factor;
  Mat pMat;
};

#endif
