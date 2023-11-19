#ifndef __INS_ENSTROPY_3D_H
#define __INS_ENSTROPY_3D_H

#include <vector>

#include "measurement_3d.h"

class Enstropy3D : public Measurement3D {
public:
  Enstropy3D(INSSolverBase3D *i, const DG_FP refMu, const DG_FP refRho, 
             const DG_FP refVel,const DG_FP volume, const int sample_iter = 20);

  virtual void measure() override;

protected:
  virtual std::string get_filename() override;
  virtual std::string get_csv_header() override;
  virtual std::string get_next_csv_line() override;

private:
  DG_FP calc_enstropy();
  DG_FP calc_ke();

  struct EnstropyHistory {
    DG_FP time;
    DG_FP enstropy;
    DG_FP ke;
  };
  std::vector<EnstropyHistory> history;
  DG_FP mu, vol, rho, vel;
};

#endif