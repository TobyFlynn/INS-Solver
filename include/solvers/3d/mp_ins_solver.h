#ifndef __INS_MP_INS_SOLVER_3D_H
#define __INS_MP_INS_SOLVER_3D_H

#include "ins_solver_base.h"

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"
#include "solvers/3d/ls_solver.h"
#include "dg_matrices/3d/factor_poisson_coarse_matrix_3d.h"
#include "dg_matrices/3d/factor_poisson_matrix_free_diag_3d.h"
#include "dg_matrices/3d/factor_poisson_matrix_free_diag_oi_3d.h"
#include "dg_matrices/3d/factor_poisson_matrix_free_block_diag_3d.h"
#include "dg_linear_solvers/linear_solver.h"
#include "dg_linear_solvers/petsc_pmultigrid.h"

#include "slip_matrix/3d/factor_viscous_matrix.h"
#include "slip_matrix/3d/viscous_solver.h"

#include <string>
#include <vector>

class MPINSSolver3D : public INSSolverBase3D {
public:
  MPINSSolver3D(DGMesh3D *m, const DG_FP re);
  MPINSSolver3D(DGMesh3D *m, const DG_FP re, const std::string &filename, const int iter);
  ~MPINSSolver3D() override;

  void step() override;
  void init() override;
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

  void apply_pressure_neumann_bc(op_dat divVelT);
  void apply_pressure_neumann_bc_oi(op_dat divVelT);
  void update_pressure_matrices(DGTempDat &pr_factor);
  void update_pressure_matrices_oi(DGTempDat &pr_factor, DGTempDat &pr_factor_oi, DGTempDat &pr_factor_surf_oi);
  void update_pressure_gradient(op_dat dpdx, op_dat dpdy, op_dat dpdz);
  void update_pressure_gradient_oi(op_dat dpdx, op_dat dpdy, op_dat dpdz);

  FactorPoissonCoarseMatrix3D *coarsePressureMatrix;
  PoissonMatrixFreeDiag *pressureMatrix;
  PoissonMatrix *viscosityMatrix;
  PETScPMultigrid *pressureSolver;
  LinearSolver *viscositySolver = nullptr;
  LevelSetSolver3D *lsSolver;

  FactorViscousMatrix3D *slipViscousMatrix;
  ViscousSolver3D *slipViscousSolver;

  DG_FP reynolds;
  bool resuming, surface_tension, pr_over_int, uses_slip_bcs;

  op_dat tmp_bc_1, tmp_npf_bc;
  op_dat pr_bc, pr_bc_types;
  op_dat vis_bc_types, vis_bc;
  op_dat art_vis, st[2][3];
  op_dat bc_data_2, vis_bc_types_2, bc_data_3, vis_bc_types_3;
};

#endif
