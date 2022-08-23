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
  void set_linear_solver(int ls);

  double getAvgPressureConvergance();
  double getAvgViscosityConvergance();

  DGMesh *mesh;
  INSData *data;
  LS *ls;
  double dt;
private:
  PetscPressureSolve *pressurePoisson;
  PetscViscositySolve *viscosityPoisson;
  PMultigrid *pMultigrid;
  int problem;
  int linear_solver;
};

#endif
