#ifndef __INS_TEMPERATURE_SOLVER_H
#define __INS_TEMPERATURE_SOLVER_H

#include "ins_solver_base.h"

#include "dg_compiler_defs.h"

#include <string>
#include <vector>

#include "dg_matrices/2d/factor_poisson_coarse_matrix_2d.h"
#include "dg_matrices/2d/factor_poisson_matrix_free_diag_2d.h"
#include "dg_matrices/2d/mm_poisson_matrix_free_2d.h"
#include "dg_linear_solvers/petsc_pmultigrid.h"
#include "dg_linear_solvers/petsc_inv_mass.h"
#include "advec_diff_solver.h"

#include "dg_mesh/dg_mesh_2d.h"

class TemperatureAdvecDiffSolver2D : public AdvecDiffSolver2D {
public:
  TemperatureAdvecDiffSolver2D(DGMesh2D *m);
  void set_bc_types(op_dat bc);

protected:
  virtual void bc_kernel_advec(op_dat val, op_dat u, op_dat v, op_dat out) override;
  virtual void bc_kernel_diff_0(op_dat val, op_dat out_x, op_dat out_y) override;
  virtual void bc_kernel_diff_1(op_dat val, op_dat val_x, op_dat val_y, op_dat vis, op_dat out) override;

  op_dat bc_types;
};

class INSTemperatureSolver2D : public INSSolverBase2D {
public:
  INSTemperatureSolver2D(DGMesh2D *m);
  INSTemperatureSolver2D(DGMesh2D *m, const std::string &filename, const int iter);
  ~INSTemperatureSolver2D();

  void init(const DG_FP re, const DG_FP refVel) override;
  void step() override;

  void dump_visualisation_data(const std::string &filename) override;
  void dump_checkpoint_data(const std::string &filename) override;
  void save_l2_err_history(const std::string &filename);

private:
  void setup_common();
  void advection();
  bool pressure();
  bool viscosity();

  void update_rho();
  void update_temperature();

  FactorPoissonCoarseMatrix2D *pressureCoarseMatrix;
  FactorPoissonMatrixFreeDiag2D *pressureMatrix;
  MMPoissonMatrixFree2D *viscosityMatrix;
  PETScPMultigrid *pressureSolver;
  PETScInvMassSolver *viscositySolver;
  TemperatureAdvecDiffSolver2D *advecDiffSolver;

  bool resuming;
  bool dt_forced;
  DG_FP reynolds;

  op_dat pr_bc_types, vis_bc_types;
  op_dat temperature, rho;
};

#endif
