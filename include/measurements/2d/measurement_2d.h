#ifndef __INS_MEASUREMENT_2D_H
#define __INS_MEASUREMENT_2D_H

#include <string>

#include "solvers/2d/ins_solver_base.h"

class Measurement2D {
public:
  Measurement2D(INSSolverBase2D *i, const int sample_iter);

  virtual void measure() = 0;
  virtual void output(std::string &filename) = 0;

protected:
  std::string double_to_text(const double &d);
  bool sample_this_iter();

  INSSolverBase2D *ins;

private:
  int sample_rate;
  int count;
};

#endif