#ifndef __INS_LIFT_DRAG_CYLINDER_2D_H
#define __INS_LIFT_DRAG_CYLINDER_2D_H

#include <vector>

#include "measurement_2d.h"
#include "solvers/2d/ins_solver_base.h"

class LiftDragCylinder2D : public Measurement2D {
public:
  LiftDragCylinder2D(INSSolverBase2D *i, const DG_FP refMu, const DG_FP x0, 
                     const DG_FP y0, const DG_FP x1, const DG_FP y1,
                     const int sample_iter = 20);

  virtual void measure() override;

protected:
  virtual std::string get_filename() override;
  virtual std::string get_csv_header() override;
  virtual std::string get_next_csv_line() override;

private:
  struct LiftDragHistory {
    DG_FP time;
    DG_FP lift;
    DG_FP drag;
  };
  std::vector<LiftDragHistory> history;
  DG_FP mu;
  DG_FP box[4];
};

#endif