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

class INSSolver3D : public INSSolverBase3D {
public:
  INSSolver3D(DGMesh3D *m);
  INSSolver3D(DGMesh3D *m, const std::string &filename, const int iter);
  ~INSSolver3D();

  void init(const DG_FP re, const DG_FP refVel) override;
  void step() override;

private:
  void setup_common();
  void advection();
  void pressure();
  void viscosity();

  // DGMesh3D *mesh;
  // PoissonMatrix3D *pressureMatrix;
  PoissonCoarseMatrix3D *pressureCoarseMatrix;
  // PoissonSemiMatrixFree3D *pressureMatrix;
  PoissonMatrixFreeDiag3D *pressureMatrix;
  // MMPoissonMatrix3D *viscosityMatrix;
  MMPoissonMatrixFree3D *viscosityMatrix;
  LinearSolver *pressureSolver;
  PETScInvMassSolver *viscositySolver;
  DG_FP reynolds;
  bool resuming;
  bool forced_dt;

  op_dat tmp_bc_1, tmp_npf_bc;
  op_dat pr_bc, pr_bc_types;
  op_dat vis_bc_types, vis_bc;
};

#endif
