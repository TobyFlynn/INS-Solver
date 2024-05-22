#ifndef __INS_ENSTROPHY_3D_H
#define __INS_ENSTROPHY_3D_H

#include <vector>

#include "measurement_3d.h"
#include "solvers/3d/ins_solver_base.h"

class Enstrophy3D : public Measurement3D {
public:
  Enstrophy3D(SimulationDriver *d, const DG_FP refMu, const DG_FP refRho,
              const DG_FP refVel, const DG_FP refLen, const DG_FP volume,
              const int sample_iter = 20);

  virtual void measure() override;

protected:
  virtual std::string get_filename() override;
  virtual std::string get_csv_header() override;
  virtual std::string get_next_csv_line() override;

private:
  DG_FP calc_enstrophy();
  DG_FP calc_ke();

  struct EnstrophyHistory {
    DG_FP time;
    DG_FP enstrophy;
    DG_FP ke;
  };
  std::vector<EnstrophyHistory> history;
  DG_FP mu, vol, rho, len, vel;

  INSSolverBase3D *ins;
};

#endif
