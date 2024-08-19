#ifndef __INS_ZALESAK_ERROR_2D_H
#define __INS_ZALESAK_ERROR_2D_H

#include <vector>

#include "measurement_2d.h"
#include "drivers/2d/ls_driver.h"

class ZalesakError : public Measurement2D {
public:
  ZalesakError(SimulationDriver *d, const DG_FP x0, const DG_FP y0, 
               const DG_FP x1, const DG_FP y1, const int sample_iter = 1);

  virtual void measure() override;

protected:
  virtual std::string get_filename() override;
  virtual std::string get_csv_header() override;
  virtual std::string get_next_csv_line() override;

private:
  struct ErrorHistory {
    DG_FP time;
    DG_FP l1;
    DG_FP l2;
    DG_FP l_max;
  };
  std::vector<ErrorHistory> history;
  LSDriver2D *ls_driver;
  DG_FP box[4];
};

#endif