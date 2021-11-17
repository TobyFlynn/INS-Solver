#ifndef __INS_SOLVER_H
#define __INS_SOLVER_H

#include <string>

#include "ins_data.h"
#include "poisson.h"
#include "ls.h"

#include "dg_mesh.h"

class Solver {
public:
  Solver(std::string filename, bool pre, int prob);
  ~Solver();

  void advection(int currentInd, double a0, double a1, double b0, double b1,
                 double g0, double t);

  bool pressure(int currentInd, double a0, double a1, double b0, double b1,
                double g0, double t);

  bool viscosity(int currentInd, double a0, double a1, double b0, double b1,
                 double g0, double t);

  void update_surface(int currentInd);

  double getAvgPressureConvergance();
  double getAvgViscosityConvergance();

  INSData *data;
  DGMesh *mesh;
  LS *ls;
  double dt;
private:
  Poisson_MF2 *pressurePoisson;
  Poisson_MF2 *viscosityPoisson;

  int problem;
};

#endif
