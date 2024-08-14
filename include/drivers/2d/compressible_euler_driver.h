#ifndef __INS_COMPRESSIBLE_EULER_DRIVER_2D_H
#define __INS_COMPRESSIBLE_EULER_DRIVER_2D_H

#include "simulation_driver.h"

#include "solvers/2d/compressible_euler.h"

class CompressibleEulerDriver2D : public SimulationDriver {
public:
  CompressibleEulerDriver2D(DGMesh2D *m);
  ~CompressibleEulerDriver2D() override;

  void init() override;
  void step() override;
  void dump_visualisation_data(const std::string &filename) override;
  void dump_checkpoint_data(const std::string &filename) override;
  DG_FP get_time() override;
  DGMesh2D* get_mesh();
  op_dat get_rho();

private:
  DGMesh2D *mesh;
  DG_FP time, h, dt;
  CompressibleEuler2D *euler_solver;
};

#endif