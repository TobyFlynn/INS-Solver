#ifndef __INS_L2_EULER_VORTEX_2D_H
#define __INS_L2_EULER_VORTEX_2D_H

#include <vector>

#include "measurement_2d.h"
#include "drivers/2d/compressible_euler_driver.h"

class L2EulerVortex2D : public Measurement2D {
public:
  L2EulerVortex2D(SimulationDriver *d, const int sample_iter = 1);

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
  CompressibleEulerDriver2D *euler_driver;
};

#endif