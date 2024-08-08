#ifndef __INS_3D_INS_SOLVER_H
#define __INS_3D_INS_SOLVER_H

#include "ins_solver_base.h"

#include "dg_compiler_defs.h"

#include "dg_mesh/dg_mesh_3d.h"

#include "dg_matrices/3d/poisson_matrix_3d.h"
#include "dg_matrices/3d/poisson_coarse_matrix_3d.h"
#include "dg_matrices/3d/poisson_matrix_free_diag_3d.h"
#include "dg_matrices/3d/poisson_matrix_free_block_diag_3d.h"
#include "dg_matrices/3d/mm_poisson_matrix_3d.h"
#include "dg_matrices/3d/mm_poisson_matrix_free_3d.h"
#include "dg_linear_solvers/petsc_inv_mass.h"
#include "dg_linear_solvers/petsc_pmultigrid.h"

#include "slip_matrix/3d/matrix_3_vec.h"
#include "slip_matrix/3d/viscous_solver.h"

class INSSolver3D : public INSSolverBase3D {
public:
  INSSolver3D(DGMesh3D *m, const DG_FP re);
  INSSolver3D(DGMesh3D *m, const DG_FP re, const std::string &filename, const int iter);
  ~INSSolver3D() override;

  void init() override;
  void step() override;

private:
  void setup_common();
  void advection();
  void pressure();
  void viscosity();

  PoissonCoarseMatrix3D *pressureCoarseMatrix;
  PoissonMatrixFreeDiag3D *pressureMatrix;
  PoissonMatrix *viscosityMatrix;
  PETScPMultigrid *pressureSolver;
  LinearSolver *viscositySolver = nullptr;
  DG_FP reynolds, fsv_relaxation_factor, fsv_factor;
  bool resuming, forced_dt, force_superficial_velocity, uses_slip_bcs;

  Matrix3Vec *slipViscousMatrix;
  ViscousSolver3D *slipViscousSolver;

  op_dat tmp_bc_1, tmp_npf_bc;
  op_dat pr_bc, pr_bc_types;
  op_dat vis_bc_types, vis_bc;
  op_dat bc_data_2, vis_bc_types_2, bc_data_3, vis_bc_types_3;
};

#endif
