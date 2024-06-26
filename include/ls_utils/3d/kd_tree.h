#ifndef __INS_KD_TREE_3D_H
#define __INS_KD_TREE_3D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include <vector>
#include <set>
#include <map>

#include "dg_mesh/dg_mesh_3d.h"
#include "poly_approx.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

struct KDCoord {
  DG_FP x_rot;
  DG_FP y_rot;
  DG_FP z_rot;
  DG_FP x;
  DG_FP y;
  DG_FP z;
  int poly;
  int rank;
};

struct KDNode {
  int l;
  int r;
  std::vector<KDCoord>::iterator start;
  std::vector<KDCoord>::iterator end;
  arma::mat rot;
  DG_FP x_min, x_max;
  DG_FP y_min, y_max;
  DG_FP z_min, z_max;
  int axis;
};

class KDTree3D {
public:
  KDTree3D(DGMesh3D *m);

  virtual void build_tree(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                          const int num, op_dat s);
  virtual void reset();
  KDCoord closest_point(DG_FP x, DG_FP y, DG_FP z);
  std::vector<PolyApprox3D> get_polys();
  void set_poly_data(std::vector<PolyApprox3D> &_polys, std::map<int,int> &_cell2polyMap);

  std::vector<KDCoord> points;
  bool empty;
protected:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, bool has_transformed, int level);
  DG_FP bb_sqr_dist(const int node_ind, const DG_FP x, const DG_FP y, const DG_FP z);
  void nearest_neighbour(DG_FP x, DG_FP y, DG_FP z, int current_ind, std::vector<KDCoord>::iterator &closest_pt, DG_FP &closest_distance);

  void update_poly_inds(std::vector<KDCoord> &points);
  void pre_build_setup(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                       const int num, op_dat s);

  std::vector<KDNode> nodes;

  int n;
  std::map<int,int> cell2polyMap;
  std::vector<PolyApprox3D> polys;
  static const int leaf_size = 100;
  DGMesh3D *mesh;
};

#endif
