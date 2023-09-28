#ifndef __INS_MP_INS_SOLVER_2D_H
#define __INS_MP_INS_SOLVER_2D_H

#include "ins_solver_base.h"

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_2d.h"
#include "solvers/2d/ls_solver.h"
#include "dg_matrices/2d/factor_poisson_coarse_matrix_2d.h"
#include "dg_matrices/2d/factor_poisson_matrix_free_diag_2d.h"
#include "dg_matrices/2d/factor_mm_poisson_matrix_free_diag_2d.h"
#include "dg_linear_solvers/petsc_pmultigrid.h"
#include "dg_linear_solvers/petsc_jacobi.h"

#include <string>
#include <vector>

class MPINSSolver2D : public INSSolverBase2D {
public:
  MPINSSolver2D(DGMesh2D *m);
  MPINSSolver2D(DGMesh2D *m, const std::string &filename, const int iter);
  ~MPINSSolver2D();

  void step() override;
  void init(const DG_FP re, const DG_FP refVel) override;
  void dump_data(const std::string &filename);

  op_dat rho, mu, st[2][2];
private:
  void setup_common();
  void advection();
  bool pressure();
  bool viscosity();
  void surface();

  FactorPoissonCoarseMatrix2D *pressureCoarseMatrix;
  FactorPoissonMatrixFreeDiag2D *pressureMatrix;
  FactorMMPoissonMatrixFreeDiag2D *viscosityMatrix;
  PETScPMultigrid *pressureSolver;
  PETScJacobiSolver *viscositySolver;
  LevelSetSolver2D *lsSolver;
  DG_FP reynolds;
  bool resuming, dt_forced, surface_tension;

  op_dat tmp_bc_1, tmp_npf_bc;
  op_dat pr_bc, pr_bc_types;
  op_dat vis_bc_types, vis_bc;
};

#endif
