#ifndef __PETSC_PMULTIGRID_H
#define __PETSC_PMULTIGRID_H

#include "op_seq.h"
#include "linear_solver.h"
#include "petscvec.h"
#include "petscksp.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "pmultigrid.h"

class PETScPMultigrid : public LinearSolver {
public:
  PETScPMultigrid(DGMesh2D *m);
  ~PETScPMultigrid();

  bool solve(op_dat rhs, op_dat ans) override;

  void calc_rhs(const double *in_d, double *out_d);
  void precond(const double *in_d, double *out_d);

private:
  void create_shell_mat();
  void set_shell_pc(PC pc);

  DGMesh2D *mesh;
  KSP ksp;

  op_dat in, out;
  bool pMatInit;
  Mat pMat;
  PMultigridPoissonSolver *pmultigridSolver;
};

#endif
