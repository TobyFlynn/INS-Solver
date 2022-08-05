#ifndef __INS_KD_TREE_H
#define __INS_KD_TREE_H

#include "op_seq.h"

#include <vector>
#include <set>
#include <map>

#include "dg_mesh.h"
#include "ls_reinit_poly.h"

struct KDCoord {
  double x;
  double y;
  int poly;
  int rank;
};

struct KDNode {
  int l;
  int r;
  KDCoord coord;
};

class KDTree {
public:
  KDTree(const double *x, const double *y, const int num, DGMesh *mesh, op_dat s);

  KDCoord closest_point(const double x, const double y);
  std::vector<PolyApprox> get_polys();

private:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, int axis);
  void nearest_neighbour(const double x, const double y, int current_ind, int axis, int &closest_ind, double &closest_distance, KDCoord &pt);

  std::set<int> cell_inds(std::vector<KDCoord> &points);
  void construct_polys(std::vector<KDCoord> &points, DGMesh *mesh, op_dat s);
  void update_poly_inds(std::vector<KDCoord> &points);

  std::vector<KDNode> nodes;
  int n;
  std::map<int,int> cell2polyMap;
  std::vector<PolyApprox> polys;
};

#endif
