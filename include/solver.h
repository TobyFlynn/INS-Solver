#ifndef __INS_SOLVER_H
#define __INS_SOLVER_H

#include <string>

#include "ins_data.h"
#include "poisson.h"

class Solver {
public:
  Solver(std::string filename, int pmethod);
  ~Solver();

  void advection(int currentInd, double a0, double a1, double b0, double b1,
                 double g0, double t);

  bool pressure(int currentInd, double a0, double a1, double b0, double b1,
                double g0, double t);

  bool viscosity(int currentInd, double a0, double a1, double b0, double b1,
                 double g0, double t);

  void lift_drag_coeff(double *lift, double *drag, int ind);

  INSData *data;
  double dt;
private:
  CubatureData *cubatureData;
  GaussData *gaussData;
  Poisson *pressurePoisson;
  Poisson *viscosityPoisson;
};

#endif
