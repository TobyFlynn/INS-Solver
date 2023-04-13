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

  op_dat tmp_np[12], tmp_npf[3], tmp_bc_1, tmp_npf_bc;
  op_dat f[3][3], n[2][3], advec_flux[3], curlVel[3];
  op_dat curl2Vel[3], dPdN[2], pr_bc, pr_bc_types, dpdx, dpdy, dpdz;
  op_dat vis_bc_types, vis_bc, bc_types;
  op_dat art_vis, vis_coeff, vis_mm, divVelT, visRHS[3];
  op_dat projRHS[3], proj_h, proj_pen;
  op_dat advec_sc[3], advec_sc_rk[3][3], advec_sc_tmp[3];
};

#endif
