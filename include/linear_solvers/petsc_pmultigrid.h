#ifndef __PETSC_PMULTIGRID_H
#define __PETSC_PMULTIGRID_H

#include "op_seq.h"
#include "linear_solver.h"
#include "petscvec.h"
#include "petscksp.h"
#include "dg_mesh/dg_mesh.h"
#include "pmultigrid.h"

class PETScPMultigrid : public LinearSolver {
public:
  PETScPMultigrid(DGMesh *m);
  ~PETScPMultigrid();

  void init() override;

  void set_coarse_matrix(PoissonCoarseMatrix *c_mat);
  bool solve(op_dat rhs, op_dat ans) override;

  void calc_rhs(const DG_FP *in_d, DG_FP *out_d);
  void precond(const DG_FP *in_d, DG_FP *out_d);

private:
  void create_shell_mat();
  void set_shell_pc(PC pc);

  DGMesh *mesh;
  KSP ksp;
  Vec b, x;

  bool pMatInit;
  Mat pMat;
  PMultigridPoissonSolver *pmultigridSolver;
};

#endif
