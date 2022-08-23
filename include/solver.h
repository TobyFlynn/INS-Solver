#ifndef __INS_SOLVER_H
#define __INS_SOLVER_H

#include <string>

#include "ins_data.h"
#include "petsc_poisson.h"
#include "p_multigrid.h"
#include "ls.h"

#include "dg_mesh.h"

class Solver {
public:
  Solver(std::string filename, int prob);
  ~Solver();

  void reverse_vel();

  void advection(int currentInd, double a0, double a1, double b0, double b1,
                 double g0, double t);

  bool pressure(int currentInd, double a0, double a1, double b0, double b1,
                double g0, double t);

  bool viscosity(int currentInd, double a0, double a1, double b0, double b1,
                 double g0, double t);

  void update_surface(int currentInd);

  double getAvgPressureConvergance();
  double getAvgViscosityConvergance();

  DGMesh *mesh;
  INSData *data;
  LS *ls;
  double dt;
private:
  void advection_non_linear(op_dat u, op_dat v, op_dat Nx, op_dat Ny, double t);
  void advection_non_linear(op_dat u0, op_dat v0, op_dat u1, op_dat v1, op_dat Nx, op_dat Ny, double t);

  PetscPressureSolve *pressurePoisson;
  PetscViscositySolve *viscosityPoisson;
  PMultigrid *pMultigrid;
  int problem;
};

#endif
