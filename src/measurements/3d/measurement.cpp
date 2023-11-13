#include "measurements/3d/measurement_3d.h"

#include <iomanip>
#include <sstream>

Measurement3D::Measurement3D(INSSolverBase3D *i, const int sample_iter) {
  ins = i;
  sample_rate = sample_iter;
  count = 0;
}

std::string Measurement3D::double_to_text(const double &d) {
    std::stringstream ss;
    ss << std::setprecision(15);
    ss << d;
    return ss.str();
}

bool Measurement3D::sample_this_iter() {
  if(count != sample_rate) {
    count++;
    return false;
  }
  count = 0;
  return true;
}