#ifndef __INS_L2_ERROR_VORTEX_2D_H
#define __INS_L2_ERROR_VORTEX_2D_H

#include <vector>

#include "measurement_2d.h"
#include "solvers/2d/ins_solver_base.h"

class L2ErrorVortex2D : public Measurement2D {
public:
  L2ErrorVortex2D(SimulationDriver *d, const int sample_iter = 5);

  virtual void measure() override;

protected:
  virtual std::string get_filename() override;
  virtual std::string get_csv_header() override;
  virtual std::string get_next_csv_line() override;

private:
  struct L2ErrorHistory {
    DG_FP time;
    DG_FP err;
  };
  std::vector<L2ErrorHistory> history;
  INSSolverBase2D *ins;
};

#endif