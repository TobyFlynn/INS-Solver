#ifndef __INS_KD_TREE_H
#define __INS_KD_TREE_H

#include "op_seq.h"

#include <vector>
#include <set>
#include <map>

#include "dg_mesh.h"
#include "ls_reinit_poly.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

struct KDCoord {
  double x_rot;
  double y_rot;
  double x;
  double y;
  int poly;
  int rank;
};

struct KDNode {
  int l;
  int r;
  std::vector<KDCoord>::iterator start;
  std::vector<KDCoord>::iterator end;
  arma::mat rot;
  double x_min, x_max;
  double y_min, y_max;
  int axis;
};

class KDTree {
public:
  KDTree(const double *x, const double *y, const int num, DGMesh *mesh, op_dat s);

  KDCoord closest_point(double x, double y);
  std::vector<PolyApprox> get_polys();

private:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, bool has_transformed, int level);
  double bb_sqr_dist(const int node_ind, const double x, const double y);
  void nearest_neighbour(double x, double y, int current_ind, std::vector<KDCoord>::iterator &closest_pt, double &closest_distance);

  std::set<int> cell_inds(std::vector<KDCoord> &points);
  void construct_polys(std::vector<KDCoord> &points, DGMesh *mesh, op_dat s);
  void update_poly_inds(std::vector<KDCoord> &points);

  std::vector<KDNode> nodes;
  std::vector<KDCoord> points;
  int n;
  std::map<int,int> cell2polyMap;
  std::vector<PolyApprox> polys;
  static const int leaf_size = 10;
};

#endif
