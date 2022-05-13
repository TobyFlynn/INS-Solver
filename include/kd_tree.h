#ifndef __INS_KD_TREE_H
#define __INS_KD_TREE_H

#include <vector>

struct KDCoord {
  double x;
  double y;
};

struct KDNode {
  int l;
  int r;
  KDCoord coord;
};

class KDTree {
public:
  KDTree(const double *x, const double *y, const int num);

  double closest_point(const double x, const double y);

private:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, int axis);
  void nearest_neighbour(const double x, const double y, int current_ind, int axis, int &closest_ind, double &closest_distance);

  std::vector<KDNode> nodes;
  int n;
};

#endif
