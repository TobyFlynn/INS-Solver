#ifndef __PETSC_BLOCK_JACOBI_H
#define __PETSC_BLOCK_JACOBI_H

#include "op_seq.h"
#include "linear_solver.h"
#include "petscvec.h"
#include "petscksp.h"
#include "dg_mesh/dg_mesh.h"

class PETScBlockJacobiSolver : public LinearSolver {
public:
  PETScBlockJacobiSolver(DGMesh *m);
  ~PETScBlockJacobiSolver();

  bool solve(op_dat rhs, op_dat ans) override;

  void calc_rhs(const DG_FP *in_d, DG_FP *out_d);
  void precond(const DG_FP *in_d, DG_FP *out_d);

private:
  void calc_precond_mat();
  void create_shell_mat();
  void set_shell_pc(PC pc);

  DGMesh *mesh;
  KSP ksp;
  op_dat pre;
  bool pMatInit;
  Mat pMat;
};

#endif
