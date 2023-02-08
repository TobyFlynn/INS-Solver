#ifndef __INS_P_MULTIGRID_H
#define __INS_P_MULTIGRID_H

#include "op_seq.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "petscvec.h"
#include "petscksp.h"
#include "linear_solver.h"
#include "petsc_amg.h"

class PMultigridPoissonSolver : public LinearSolver {
public:
  PMultigridPoissonSolver(DGMesh2D *m);
  ~PMultigridPoissonSolver();

  bool solve(op_dat rhs, op_dat ans) override;

  void calc_rhs(const double *u_d, double *rhs_d);
private:
  void cycle(int order);

  double maxEigenValue();
  void setRandomVector(op_dat vec);
  void setupDirectSolve();

  DGMesh2D *mesh;
  PETScAMGSolver *coarseSolver;

  op_dat tmp_dat[DG_ORDER], u_dat[DG_ORDER], b_dat[DG_ORDER];
};

#endif
