#ifndef __INS_MEASUREMENT_2D_H
#define __INS_MEASUREMENT_2D_H

#include <string>

class Measurement2D {
public:
  virtual void measure() = 0;
  virtual void output(std::string &filename) = 0;
};

#endif