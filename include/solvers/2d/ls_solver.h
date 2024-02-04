#ifndef __INS_LS_2D_H
#define __INS_LS_2D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include "solvers/2d/advection_solver.h"
#include "ls_utils/2d/ls_reinit_poly.h"
#include "ls_utils/2d/kd_tree.h"

#include "dg_mesh/dg_mesh_2d.h"

#include <map>
#include <vector>

class LevelSetAdvectionSolver2D : public AdvectionSolver2D {
public:
  LevelSetAdvectionSolver2D(DGMesh2D *m);
  void set_bc_types(op_dat bc);

protected:
  virtual void bc_kernel(op_dat val, op_dat u, op_dat v, op_dat out) override;
  virtual void bc_kernel_oi(op_dat val, op_dat u, op_dat v, op_dat flux) override;

  op_dat bc_types;
};

class LevelSetSolver2D {
public:
  LevelSetSolver2D(DGMesh2D *m);
  LevelSetSolver2D(DGMesh2D *m, const std::string &filename);
  ~LevelSetSolver2D();

  void init();
  void set_bc_types(op_dat bc);

  void setVelField(op_dat u1, op_dat v1);
  void step(const DG_FP dt, const int num_steps);
  void getRhoMu(op_dat rho, op_dat mu);
  void getRhoVolOI(op_dat rho);
  void getMuVolOI(op_dat mu);
  void getRhoSurfOI(op_dat rho);
  void getNormalsCurvature(op_dat nx, op_dat ny, op_dat curv);

  DGMesh2D *mesh;

  op_dat u, v, s, dsdx, dsdy, s_sample_x, s_sample_y;

  DG_FP alpha, order_width, ls_cap;
private:
  void sampleInterface(op_dat sampleX, op_dat sampleY, std::vector<PolyApprox> &polys, 
                       std::map<int,int> &cell2polyMap, std::set<int> &cellInds);
  void reinitLS();
  bool reinitNeeded();

  DG_FP h, epsilon, reinit_dt, reinit_width;
  int numSteps;
  bool resuming;
  int reinit_counter, reinit_frequency;

  LevelSetAdvectionSolver2D *advecSolver;
  KDTree *kdtree;
};

#endif
