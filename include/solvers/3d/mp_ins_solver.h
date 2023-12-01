#ifndef __INS_MP_INS_SOLVER_3D_H
#define __INS_MP_INS_SOLVER_3D_H

#include "ins_solver_base.h"

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"
#include "solvers/3d/ls_solver.h"
#include "dg_matrices/3d/factor_poisson_matrix_3d.h"
#include "dg_matrices/3d/factor_poisson_coarse_matrix_3d.h"
#include "dg_matrices/3d/factor_poisson_semi_matrix_free_3d.h"
#include "dg_matrices/3d/factor_poisson_matrix_free_diag_3d.h"
#include "dg_matrices/3d/factor_poisson_matrix_free_diag_oi_3d.h"
#include "dg_matrices/3d/factor_mm_poisson_matrix_3d.h"
#include "dg_matrices/3d/factor_mm_poisson_matrix_free_3d.h"
#include "dg_matrices/3d/factor_mm_poisson_semi_matrix_free_3d.h"
#include "dg_matrices/3d/factor_mm_poisson_matrix_free_diag_3d.h"
#include "dg_linear_solvers/linear_solver.h"
#include "dg_linear_solvers/petsc_pmultigrid.h"
#include "dg_linear_solvers/petsc_inv_mass.h"
#include "dg_linear_solvers/petsc_jacobi.h"

#include <string>
#include <vector>

class MPINSSolver3D : public INSSolverBase3D {
public:
  MPINSSolver3D(DGMesh3D *m);
  MPINSSolver3D(DGMesh3D *m, const std::string &filename, const int iter);
  ~MPINSSolver3D();

  void step() override;
  void init(const DG_FP re, const DG_FP refVel) override;
  void dump_visualisation_data(const std::string &filename) override;
  void dump_checkpoint_data(const std::string &filename) override;

  op_dat rho, mu;
private:
  void setup_common();
  void advection();
  void pressure();
  void viscosity();
  void surface();

  void surface_tension_grad(op_dat dx, op_dat dy, op_dat dz);
  void surface_tension_curvature(op_dat curv);

  FactorPoissonCoarseMatrix3D *coarsePressureMatrix;
  // FactorPoissonSemiMatrixFree3D *pressureMatrix;
  // FactorPoissonMatrixFreeDiag3D *pressureMatrix;
  FactorPoissonMatrixFreeDiagOI3D *pressureMatrix;
  // FactorMMPoissonMatrix3D *viscosityMatrix;
  // FactorMMPoissonSemiMatrixFree3D *viscosityMatrix;
  FactorMMPoissonMatrixFreeDiag3D *viscosityMatrix;
  // LinearSolver *pressureSolver;
  PETScPMultigrid *pressureSolver;
  // LinearSolver *viscositySolver;
  PETScJacobiSolver *viscositySolver;
  LevelSetSolver3D *lsSolver;
  DG_FP reynolds;
  bool resuming, surface_tension;

  op_dat tmp_bc_1, tmp_npf_bc;
  op_dat pr_bc, pr_bc_types;
  op_dat vis_bc_types, vis_bc;
  op_dat art_vis, st[2][3];
};

#endif
