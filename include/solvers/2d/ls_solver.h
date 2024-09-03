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
  void sampleInterface(op_dat sampleX, op_dat sampleY, std::vector<PolyApprox> &polys,
                       std::map<int,int> &cell2polyMap, std::set<int> &cellInds);

  DGMesh2D *mesh;

  op_dat u, v, s, dsdx, dsdy, s_sample_x, s_sample_y, kink;

  DG_FP alpha, order_width, ls_cap, h;
private:
  void reinitLS();
  bool reinitNeeded();
  void detect_kinks();
  void create_point_map_for_kink_detection();

  DG_FP epsilon, reinit_dt, reinit_width;
  DG_FP kink_max_distance_between_points, kink_sqr_tol;
  int kink_max_neighbours;
  int numSteps;
  bool resuming, reinitialise, kink_detection, kink_avoid_whole_element;
  int reinit_counter, reinit_frequency;

  LevelSetAdvectionSolver2D *advecSolver;
  KDTree *kdtree;
  std::map<DGUtils::Vec<2>,std::vector<DGUtils::Vec<2>>> point_map_for_kink_detection;
};

#endif
