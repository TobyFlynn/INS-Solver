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

private:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, int axis);

  std::vector<KDNode> nodes;
  int n;
};

#endif
