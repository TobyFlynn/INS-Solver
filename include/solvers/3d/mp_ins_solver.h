#ifndef __INS_MP_INS_SOLVER_3D_H
#define __INS_MP_INS_SOLVER_3D_H

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"
#include "solvers/3d/ls_solver.h"
#include "matrices/3d/factor_poisson_matrix_3d.h"
#include "matrices/3d/factor_poisson_coarse_matrix_3d.h"
#include "matrices/3d/factor_poisson_semi_matrix_free_3d.h"
#include "matrices/3d/factor_mm_poisson_matrix_3d.h"
#include "matrices/3d/factor_mm_poisson_matrix_free_3d.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/petsc_pmultigrid.h"

#include <string>

class MPINSSolver3D {
public:
  MPINSSolver3D(DGMesh3D *m);
  ~MPINSSolver3D();

  void step();
  void init(const DG_FP re, const DG_FP refVel);
  DG_FP get_time();
  DG_FP get_dt();
  void dump_data(const std::string &filename);

  op_dat vel[2][3], velT[3], velTT[3], pr, rho, mu;
private:
  void advection();
  void pressure();
  void viscosity();
  void surface();
  void shock_capturing();

  DGMesh3D *mesh;
  FactorPoissonCoarseMatrix3D *coarsePressureMatrix;
  FactorPoissonSemiMatrixFree3D *pressureMatrix;
  FactorMMPoissonMatrix3D *viscosityMatrix;
  // LinearSolver *pressureSolver;
  PETScPMultigrid *pressureSolver;
  LinearSolver *viscositySolver;
  LevelSetSolver3D *lsSolver;
  DG_FP g0, a0, a1, b0, b1, dt, time, h;
  DG_FP reynolds;
  int currentInd;

  op_dat tmp_np[9], tmp_npf[3], tmp_bc_1, tmp_npf_bc;
  op_dat f[3][3], n[2][3], advec_flux[3], curlVel[3], divVelT;
  op_dat curl2Vel[3], dPdN[2], pr_bc, pr_bc_types, dpdx, dpdy, dpdz;
  op_dat vis_bc_types, vis_bc, bc_types, pr_factor, vis_factor, vis_mm_factor;
  op_dat art_vis, shock_u, shock_u_hat, shock_u_modal;
};

#endif
