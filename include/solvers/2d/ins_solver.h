#ifndef __INS_SOLVER_H
#define __INS_SOLVER_H

#include "dg_compiler_defs.h"

#include <string>

#include "matrices/2d/poisson_coarse_matrix_2d.h"
#include "matrices/2d/poisson_matrix_free_diag_2d.h"
// #include "matrices/2d/mm_poisson_matrix_over_int_2d.h"
// #include "matrices/2d/mm_poisson_matrix_free_over_int_2d.h"
// #include "matrices/2d/cub_poisson_matrix_2d.h"
// #include "matrices/2d/cub_mm_poisson_matrix_2d.h"
#include "linear_solvers/linear_solver.h"
// #include "linear_solvers/petsc_inv_mass.h"

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

  DGMesh2D *mesh;
private:
  void setup_common();
  void advection();
  bool pressure();
  void project_velocity();
  bool viscosity();

  PoissonCoarseMatrix2D *pressureCoarseMatrix;
  PoissonMatrixFreeDiag2D *pressureMatrix;
  // MMPoissonMatrixOverInt2D *viscosityMatrix;
  // MMPoissonMatrixFreeOverInt2D *viscosityMatrix;
  LinearSolver *pressureSolver;
  // LinearSolver *viscositySolver;
  // PETScInvMassSolver *viscositySolver;

  bool resuming;
  int currentInd;
  DG_FP a0, a1, b0, b1, g0, dt, time;
  DG_FP reynolds;

  op_dat vel[2][2], n[2][2], velT[2], velTT[2], pr, dPdN[2];
  op_dat bc_types, pr_bc_types, vis_bc_types;
};

#endif
