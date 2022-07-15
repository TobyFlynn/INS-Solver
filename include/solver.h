#ifndef __INS_SOLVER_H
#define __INS_SOLVER_H

#include <string>

#include "ins_data.h"
#include "poisson.h"
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
  PressureSolve *pressurePoisson;
  ViscositySolve *viscosityPoisson;
  int problem;
};

#endif
