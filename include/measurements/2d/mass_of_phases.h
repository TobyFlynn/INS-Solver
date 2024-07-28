#ifndef __INS_MASS_OF_PHASES_2D_H
#define __INS_MASS_OF_PHASES_2D_H

#include <vector>

#include "measurement_2d.h"
#include "solvers/2d/ins_solver_base.h"
#include "solvers/2d/mp_ins_solver.h"

class MassOfPhases2D : public Measurement2D {
public:
  MassOfPhases2D(SimulationDriver *d, const int sample_iter = 5);

  virtual void measure() override;

protected:
  virtual std::string get_filename() override;
  virtual std::string get_csv_header() override;
  virtual std::string get_next_csv_line() override;

private:
  struct MassHistory {
    DG_FP time;
    DG_FP phase0_vol;
    DG_FP phase1_vol;
    DG_FP phase0_mass;
    DG_FP phase1_mass;
    DG_FP total_mass;
  };
  std::vector<MassHistory> history;
  MPINSSolver2D *mpins;
};

#endif