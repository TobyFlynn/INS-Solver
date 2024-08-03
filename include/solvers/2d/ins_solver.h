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
#include "dg_linear_solvers/petsc_pmultigrid.h"
#include "dg_linear_solvers/petsc_inv_mass.h"

#include "slip_matrix/2d/matrix_2_vec.h"
#include "slip_matrix/2d/viscous_solver.h"

#include "dg_mesh/dg_mesh_2d.h"

class INSSolver2D : public INSSolverBase2D {
public:
  INSSolver2D(DGMesh2D *m, const DG_FP re);
  INSSolver2D(DGMesh2D *m, const DG_FP re, const std::string &filename, const int iter);
  ~INSSolver2D() override;

  void init() override;
  void step() override;

private:
  void setup_common();
  void advection();
  bool pressure();
  bool viscosity();
  void no_viscosity();

  PoissonCoarseMatrix2D *pressureCoarseMatrix;
  PoissonMatrixFreeDiag2D *pressureMatrix;
  PoissonMatrix *viscosityMatrix;
  PETScPMultigrid *pressureSolver;
  LinearSolver *viscositySolver = nullptr;

  Matrix2Vec *slipViscousMatrix;
  ViscousSolver *slipViscousSolver;

  bool resuming;
  bool vis_solve;
  bool dt_forced;
  bool uses_slip_bcs;
  DG_FP reynolds;

  op_dat pr_bc_types, vis_bc_types, bc_data_2, vis_bc_types_2;

  std::vector<std::pair<DG_FP,DG_FP>> l2_err_history;
};

#endif
