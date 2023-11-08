#ifndef __INS_SOLVER_H
#define __INS_SOLVER_H

#include "ins_solver_base.h"

#include "dg_compiler_defs.h"

#include <string>
#include <vector>

#include "dg_matrices/2d/poisson_coarse_matrix_2d.h"
#include "dg_matrices/2d/poisson_matrix_free_diag_2d.h"
#include "dg_matrices/2d/mm_poisson_matrix_free_2d.h"
#include "dg_linear_solvers/linear_solver.h"
#include "dg_linear_solvers/petsc_inv_mass.h"

#include "dg_mesh/dg_mesh_2d.h"

class INSSolver2D : public INSSolverBase2D {
public:
  INSSolver2D(DGMesh2D *m);
  INSSolver2D(DGMesh2D *m, const std::string &filename, const int iter);
  ~INSSolver2D();

  void init(const DG_FP re, const DG_FP refVel) override;
  void step() override;

  void save_l2_err_history(const std::string &filename);

private:
  void setup_common();
  void advection();
  bool pressure();
  bool viscosity();
  void no_viscosity();

  DG_FP l2_vortex_error(DG_FP time);
  void record_l2_err();

  PoissonCoarseMatrix2D *pressureCoarseMatrix;
  PoissonMatrixFreeDiag2D *pressureMatrix;
  MMPoissonMatrixFree2D *viscosityMatrix;
  LinearSolver *pressureSolver;
  PETScInvMassSolver *viscositySolver;

  bool resuming;
  bool vis_solve;
  bool dt_forced;
  DG_FP reynolds;

  op_dat pr_bc_types, vis_bc_types;

  std::vector<std::pair<DG_FP,DG_FP>> l2_err_history;
};

#endif
