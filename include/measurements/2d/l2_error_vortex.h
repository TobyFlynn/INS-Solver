#ifndef __INS_L2_ERROR_VORTEX_2D_H
#define __INS_L2_ERROR_VORTEX_2D_H

#include <vector>

#include "measurement_2d.h"
#include "solvers/2d/ins_solver_base.h"

class L2ErrorVortex2D : public Measurement2D {
public:
  L2ErrorVortex2D(INSSolverBase2D *i, const int sample_iter = 20);

  virtual void measure() override;
  virtual void output(std::string &path) override;

private:
  struct L2ErrorHistory {
    DG_FP time;
    DG_FP err;
  };
  std::vector<L2ErrorHistory> history;
};

#endif