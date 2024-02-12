#ifndef __INS_LS_DRIVER_2D_H
#define __INS_LS_DRIVER_2D_H

#include "simulation_driver.h"
#include "solvers/2d/ls_solver.h"

class LSDriver2D : public SimulationDriver {
public:
  LSDriver2D(DGMesh2D *m);
  ~LSDriver2D() override;

  void init() override;
  void step() override;
  void dump_visualisation_data(const std::string &filename) override;
  void dump_checkpoint_data(const std::string &filename) override;
  DG_FP get_time() override;
  DGMesh2D* get_mesh();
  op_dat get_surface();

private:
  DGMesh2D *mesh;
  LevelSetSolver2D *lsSolver;
  op_dat u, v, bc_types;
  DG_FP time, h, dt;
};

#endif