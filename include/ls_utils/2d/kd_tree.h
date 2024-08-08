#ifndef __INS_KD_TREE_H
#define __INS_KD_TREE_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include <vector>
#include <set>
#include <map>

#include "dg_utils.h"
#include "dg_mesh/dg_mesh_2d.h"
#include "ls_reinit_poly.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

struct KDCoord {
  DGUtils::Vec<2> coord_rot;
  DGUtils::Vec<2> coord;
  int poly;
  int rank;
};

struct KDNode {
  int l;
  int r;
  std::vector<KDCoord>::iterator start;
  std::vector<KDCoord>::iterator end;
  arma::mat rot;
  DGUtils::Vec<2> min;
  DGUtils::Vec<2> max;
  int axis;
};

class KDTree {
public:
  KDTree(DGMesh2D *m);

  virtual void build_tree(const DG_FP *x, const DG_FP *y, const int num, op_dat s);
  virtual void reset();
  KDCoord closest_point(DG_FP x, DG_FP y);
  std::vector<PolyApprox> get_polys();
  void set_poly_data(std::vector<PolyApprox> &_polys, std::map<int,int> &_cell2polyMap);

  std::vector<KDCoord> points;
  bool empty;
protected:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, bool has_transformed, int level);
  DG_FP bb_sqr_dist(const int node_ind, DGUtils::Vec<2> &coord);
  void nearest_neighbour(DGUtils::Vec<2> coord, int current_ind, std::vector<KDCoord>::iterator &closest_pt, DG_FP &closest_distance);

  void update_poly_inds(std::vector<KDCoord> &points);
  void pre_build_setup(const DG_FP *x, const DG_FP *y, const int num, op_dat s);

  std::vector<KDNode> nodes;

  int n;
  std::map<int,int> cell2polyMap;
  std::vector<PolyApprox> polys;
  static const int leaf_size = 10;
  DGMesh2D *mesh;
};

#endif
