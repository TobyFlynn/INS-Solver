#ifndef __INS_3D_INS_SOLVER_H
#define __INS_3D_INS_SOLVER_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"

#include "matrices/3d/poisson_matrix_3d.h"
#include "matrices/3d/poisson_coarse_matrix_3d.h"
#include "matrices/3d/poisson_semi_matrix_free_3d.h"
#include "matrices/3d/poisson_matrix_free_diag_3d.h"
#include "matrices/3d/mm_poisson_matrix_3d.h"
#include "matrices/3d/mm_poisson_matrix_free_3d.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/petsc_inv_mass.h"

class INSSolver3D {
public:
  INSSolver3D(DGMesh3D *m);
  INSSolver3D(DGMesh3D *m, const std::string &filename, const int iter);
  ~INSSolver3D();

  void init(const DG_FP re, const DG_FP refVel);
  void step();

  DG_FP get_time();
  DG_FP get_dt();
  void dump_data(const std::string &filename);

  op_dat vel[2][3], velT[3], velTT[3], pr;
private:
  void setup_common();
  void advection();
  void pressure();
  void viscosity();
  void project_velocity();
  void advec_current_non_linear();
  void advec_standard();
  void advec_sub_cycle();
  void advec_sub_cycle_rhs(op_dat u_in, op_dat v_in, op_dat w_in,
                           op_dat u_out, op_dat v_out, op_dat w_out,
                           const double t);
  void advec_sub_cycle_rk_step(const DG_FP time_sc, op_dat u, op_dat v, op_dat w);
  DG_FP max_vel();
  // void shock_capturing();

  DGMesh3D *mesh;
  // PoissonMatrix3D *pressureMatrix;
  PoissonCoarseMatrix3D *pressureCoarseMatrix;
  // PoissonSemiMatrixFree3D *pressureMatrix;
  PoissonMatrixFreeDiag3D *pressureMatrix;
  // MMPoissonMatrix3D *viscosityMatrix;
  MMPoissonMatrixFree3D *viscosityMatrix;
  LinearSolver *pressureSolver;
  // LinearSolver *viscositySolver;
  PETScInvMassSolver *viscositySolver;
  DG_FP g0, a0, a1, b0, b1, dt, sub_cycle_dt, time, h;
  DG_FP reynolds;
  int currentInd, sub_cycles, it_pre_sub_cycle;
  bool resuming;
  bool div_div_proj;

  op_dat tmp_bc_1, tmp_npf_bc;
  op_dat n[2][3];
  op_dat dPdN[2], pr_bc, pr_bc_types;
  op_dat vis_bc_types, vis_bc, bc_types;
  // op_dat art_vis, vis_coeff, vis_mm;
  op_dat proj_h;
};

#endif
