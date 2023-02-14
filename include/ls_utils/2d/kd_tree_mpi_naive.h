#ifndef __INS_KD_TREE_H
#define __INS_KD_TREE_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include <vector>
#include <set>
#include <map>

#include "dg_mesh.h"
#include "ls_reinit_poly.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

struct KDCoord {
  DG_FP x_rot;
  DG_FP y_rot;
  DG_FP x;
  DG_FP y;
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
  int axis;
};

class KDTreeMPINaive {
public:
  KDTreeMPINaive(const DG_FP *x, const DG_FP *y, const int num, DGMesh2D *mesh, op_dat s);

  KDCoord closest_point(DG_FP x, DG_FP y);
  std::vector<PolyApprox> get_polys();

private:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, bool has_transformed, int level);
  DG_FP bb_sqr_dist(const int node_ind, const DG_FP x, const DG_FP y);
  void nearest_neighbour(DG_FP x, DG_FP y, int current_ind, std::vector<KDCoord>::iterator &closest_pt, DG_FP &closest_distance);

  std::set<int> cell_inds(std::vector<KDCoord> &points);
  void construct_polys(std::vector<KDCoord> &points, DGMesh2D *mesh, op_dat s);
  void update_poly_inds(std::vector<KDCoord> &points);

  std::vector<KDNode> nodes;
  std::vector<KDCoord> points;
  int n;
  std::map<int,int> cell2polyMap;
  std::vector<PolyApprox> polys;
  static const int leaf_size = 10;
};

#endif
