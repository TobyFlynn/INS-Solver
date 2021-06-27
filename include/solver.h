#ifndef __INS_SOLVER_H
#define __INS_SOLVER_H

#include <string>

#include "ins_data.h"
#include "pressure_solve.h"
#include "viscosity_solve.h"
#include "ls.h"

class Solver {
public:
  Solver(std::string filename, int pmethod, int prob, bool multi);
  ~Solver();

  void advection(int currentInd, double a0, double a1, double b0, double b1,
                 double g0, double t);

  bool pressure(int currentInd, double a0, double a1, double b0, double b1,
                double g0, double t);

  bool viscosity(int currentInd, double a0, double a1, double b0, double b1,
                 double g0, double t);

  void update_surface(int currentInd);

  void lift_drag_coeff(double *lift, double *drag, int ind);

  double getAvgPressureConvergance();
  double getAvgViscosityConvergance();

  INSData *data;
  CubatureData *cubatureData;
  LS *ls;
  double dt;
private:
  // CubatureData *cubatureData;
  GaussData *gaussData;
  PressureSolve *pressurePoisson;
  ViscositySolve *viscosityPoisson;

  int problem;
  bool multiphase;
};

#endif
