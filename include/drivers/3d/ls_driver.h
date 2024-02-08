#ifndef __INS_LS_DRIVER_3D_H
#define __INS_LS_DRIVER_3D_H

#include "simulation_driver.h"
#include "solvers/3d/ls_solver.h"

class LSDriver3D : public SimulationDriver {
public:
  LSDriver3D(DGMesh3D *m);
  ~LSDriver3D() override;

  void init() override;
  void step() override;
  void dump_visualisation_data(const std::string &filename) override;
  void dump_checkpoint_data(const std::string &filename) override;

private:
  DGMesh3D *mesh;
  LevelSetSolver3D *lsSolver;
  op_dat u, v, w, bc_types;
  DG_FP time, h, dt;
};

#endif