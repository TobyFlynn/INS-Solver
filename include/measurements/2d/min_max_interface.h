#ifndef __INS_MIN_MAX_INTERFACE_2D_H
#define __INS_MIN_MAX_INTERFACE_2D_H

#include <vector>

#include "measurement_2d.h"
#include "solvers/2d/ins_solver_base.h"
#include "solvers/2d/mp_ins_solver.h"

class MinMaxInterface2D : public Measurement2D {
public:
  MinMaxInterface2D(INSSolverBase2D *i, const int sample_iter = 20);

  virtual void measure() override;

protected:
  virtual std::string get_filename() override;
  virtual std::string get_csv_header() override;
  virtual std::string get_next_csv_line() override;

private:
  struct MinMaxHistory {
    DG_FP time;
    DG_FP min_x;
    DG_FP min_y;
    DG_FP max_x;
    DG_FP max_y;
  };
  std::vector<MinMaxHistory> history;
  MPINSSolver2D *mpins;
};

#endif