#ifndef __INS_SOLVER_H
#define __INS_SOLVER_H

#include <string>

#include "matrices/2d/poisson_matrix_2d.h"
#include "matrices/2d/mm_poisson_matrix_2d.h"
#include "linear_solvers/linear_solver.h"
#include "solvers/2d/ls_solver.h"

#include "dg_mesh/dg_mesh_2d.h"

class INSSolver2D {
public:
  INSSolver2D(DGMesh2D *m);
  ~INSSolver2D();

  void init(const double re, const double refVel);
  void step();

  void dump_data(const std::string &filename);

  DGMesh2D *mesh;
private:
  void advection();
  bool pressure();
  bool viscosity();

  PoissonMatrix2D *pressureMatrix;
  MMPoissonMatrix2D *viscosityMatrix;
  LinearSolver *pressureSolver;
  LinearSolver *viscositySolver;

  int currentInd;
  double a0, a1, b0, b1, g0, dt, time;
  double reynolds;

  op_dat vel[2][2], n[2][2], velT[2], velTT[2], pr, dPdN[2];
  op_dat tmp_np[4], f[4], divVelT, curlVel, gradCurlVel[2], pRHS;
  op_dat dpdx, dpdy, visRHS[2];
  op_dat tmp_g_np[4], gVel[2], gAdvecFlux[2], gN[2], gGradCurl[2], prBC;
  op_dat visBC[2];
  op_dat bc_types, pr_bc_types, vis_bc_types;
};

#endif
