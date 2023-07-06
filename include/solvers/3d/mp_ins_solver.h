#ifndef __INS_MP_INS_SOLVER_3D_H
#define __INS_MP_INS_SOLVER_3D_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"
#include "dg_dat_pool.h"
#include "solvers/3d/ls_solver.h"
#include "matrices/3d/factor_poisson_matrix_3d.h"
#include "matrices/3d/factor_poisson_coarse_matrix_3d.h"
#include "matrices/3d/factor_poisson_semi_matrix_free_3d.h"
#include "matrices/3d/factor_poisson_matrix_free_diag_3d.h"
#include "matrices/3d/factor_mm_poisson_matrix_3d.h"
#include "matrices/3d/factor_mm_poisson_matrix_free_3d.h"
#include "matrices/3d/factor_mm_poisson_semi_matrix_free_3d.h"
#include "matrices/3d/factor_mm_poisson_matrix_free_diag_3d.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/petsc_pmultigrid.h"
#include "linear_solvers/petsc_inv_mass.h"
#include "linear_solvers/petsc_jacobi.h"

#include <string>
#include <vector>

class MPINSSolver3D {
public:
  MPINSSolver3D(DGMesh3D *m);
  MPINSSolver3D(DGMesh3D *m, const std::string &filename, const int iter);
  ~MPINSSolver3D();

  void step();
  void init(const DG_FP re, const DG_FP refVel);
  DG_FP get_time();
  DG_FP get_dt();
  void dump_data(const std::string &filename);

  op_dat vel[2][3], velT[3], velTT[3], pr, rho, mu;
private:
  void setup_common();
  void advection();
  void pressure();
  void viscosity();
  void surface();
  void shock_capture_filter_dat(op_dat in);
  void project_velocity();
  void advec_current_non_linear();
  void advec_standard();
  void advec_sub_cycle();
  void advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat w_in,
                           op_dat u_out, op_dat v_out, op_dat w_out,
                           const double t);
  void advec_sub_cycle_rk_step(const DG_FP time_sc, op_dat u, op_dat v, op_dat w);
  DG_FP max_vel();
  void add_to_pr_history();

  DGMesh3D *mesh;
  FactorPoissonCoarseMatrix3D *coarsePressureMatrix;
  // FactorPoissonSemiMatrixFree3D *pressureMatrix;
  FactorPoissonMatrixFreeDiag3D *pressureMatrix;
  // FactorMMPoissonMatrix3D *viscosityMatrix;
  // FactorMMPoissonSemiMatrixFree3D *viscosityMatrix;
  FactorMMPoissonMatrixFreeDiag3D *viscosityMatrix;
  // LinearSolver *pressureSolver;
  PETScPMultigrid *pressureSolver;
  // LinearSolver *viscositySolver;
  PETScJacobiSolver *viscositySolver;
  LevelSetSolver3D *lsSolver;
  DG_FP g0, a0, a1, b0, b1, dt, sub_cycle_dt, time, h;
  DG_FP reynolds;
  int currentInd, sub_cycles, it_pre_sub_cycle;
  bool resuming;
  bool div_div_proj;
  bool extrapolate_initial_guess;
  bool shock_cap;

  op_dat tmp_bc_1, tmp_npf_bc;
  op_dat n[2][3];
  op_dat dPdN[2], pr_bc, pr_bc_types;
  op_dat vis_bc_types, vis_bc, bc_types;
  op_dat art_vis;
  op_dat proj_h;

  std::vector<std::pair<DG_FP,DGTempDat>> pr_history;
};

#endif
