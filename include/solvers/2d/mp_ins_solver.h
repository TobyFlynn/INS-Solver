#ifndef __INS_MP_INS_SOLVER_2D_H
#define __INS_MP_INS_SOLVER_2D_H

#include "ins_solver_base.h"

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_2d.h"
#include "solvers/2d/ls_solver.h"
#include "solvers/2d/diffusion_solver.h"
#include "dg_matrices/2d/factor_poisson_coarse_matrix_2d.h"
#include "dg_matrices/2d/factor_poisson_matrix_free_diag_2d.h"
#include "dg_matrices/2d/factor_poisson_matrix_free_diag_oi_2d.h"
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

  void dump_visualisation_data(const std::string &filename) override;
  void dump_checkpoint_data(const std::string &filename) override;

  // Getters (used for measurements)
  op_dat get_ls();

  op_dat rho, mu, st[2][2];
private:
  void setup_common();
  void advection();
  bool pressure();
  bool viscosity();
  void surface();

  void calc_art_vis(op_dat in0, op_dat in1, op_dat out);
  void surface_tension_grad(op_dat dx, op_dat dy);
  void surface_tension_curvature(op_dat curv);

  FactorPoissonCoarseMatrix2D *pressureCoarseMatrix;
  FactorPoissonMatrixFreeDiagOI2D *pressureMatrix;
  // FactorPoissonMatrixFreeDiag2D *pressureMatrix;
  FactorMMPoissonMatrixFreeDiag2D *viscosityMatrix;
  PETScPMultigrid *pressureSolver;
  PETScJacobiSolver *viscositySolver;
  LevelSetSolver2D *lsSolver;
  DiffusionSolver2D *diffSolver;
  DG_FP reynolds;
  bool resuming, dt_forced, surface_tension;
  DG_FP st_max_diff, st_diff_width_fact, st_trans_width_fact;

  op_dat tmp_bc_1, tmp_npf_bc;
  op_dat pr_bc, pr_bc_types;
  op_dat vis_bc_types, vis_bc;
  op_dat nodes_data, nodes_count;
};

#endif
