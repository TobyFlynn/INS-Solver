#include "measurements/2d/measurement_2d.h"

#include <iomanip>
#include <sstream>

Measurement2D::Measurement2D(INSSolverBase2D *i, const int sample_iter) {
  ins = i;
  sample_rate = sample_iter;
  count = 0;
}

std::string Measurement2D::double_to_text(const double &d) {
    std::stringstream ss;
    ss << std::setprecision(15);
    ss << d;
    return ss.str();
}

bool Measurement2D::sample_this_iter() {
  if(count != sample_rate) {
    count++;
    return false;
  }
  count = 0;
  return true;
}