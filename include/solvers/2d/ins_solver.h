#ifndef __INS_SOLVER_H
#define __INS_SOLVER_H

#include "dg_compiler_defs.h"

#include <string>
#include <vector>

#include "matrices/2d/poisson_coarse_matrix_2d.h"
#include "matrices/2d/poisson_matrix_free_diag_2d.h"
#include "matrices/2d/mm_poisson_matrix_free_2d.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/petsc_inv_mass.h"

#include "dg_mesh/dg_mesh_2d.h"

class INSSolver2D {
public:
  INSSolver2D(DGMesh2D *m);
  INSSolver2D(DGMesh2D *m, const std::string &filename, const int iter);
  ~INSSolver2D();

  void init(const DG_FP re, const DG_FP refVel);
  void step();

  DG_FP get_time();
  DG_FP get_dt();
  void dump_data(const std::string &filename);
  void save_l2_err_history(const std::string &filename);

  DGMesh2D *mesh;
private:
  void setup_common();
  void advection();
  bool pressure();
  void project_velocity();
  bool viscosity();
  void no_viscosity();

  DG_FP max_vel();
  DG_FP l2_vortex_error(DG_FP time);
  void record_l2_err();

  void advec_current_non_linear();
  void advec_standard();
  void advec_sub_cycle();
  void advec_sub_cycle_rk_step(const DG_FP time_sc, op_dat u, op_dat v);
  void advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat u_out, op_dat v_out, const double t);

  PoissonCoarseMatrix2D *pressureCoarseMatrix;
  PoissonMatrixFreeDiag2D *pressureMatrix;
  MMPoissonMatrixFree2D *viscosityMatrix;
  LinearSolver *pressureSolver;
  PETScInvMassSolver *viscositySolver;

  bool resuming;
  bool div_div_proj;
  bool vis_solve;
  bool dt_forced;
  int currentInd, sub_cycles, it_pre_sub_cycle;
  DG_FP a0, a1, b0, b1, g0, dt, time, sub_cycle_dt, h;
  DG_FP reynolds;

  op_dat vel[2][2], n[2][2], velT[2], velTT[2], pr, dPdN[2];
  op_dat bc_types, pr_bc_types, vis_bc_types;
  op_dat proj_h;

  std::vector<std::pair<DG_FP,DG_FP>> l2_err_history;
};

#endif
