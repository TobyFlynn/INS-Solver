#ifndef __INS_SIMULATION_DRIVER_H
#define __INS_SIMULATION_DRIVER_H

#include <string>

class SimulationDriver {
public:
  virtual ~SimulationDriver() = default;
  virtual void init() = 0;
  virtual void step() = 0;
  virtual void dump_visualisation_data(const std::string &filename) = 0;
  virtual void dump_checkpoint_data(const std::string &filename) = 0;
};

#endif