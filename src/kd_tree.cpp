#include "kd_tree.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

bool compareX(KDCoord a, KDCoord b) {
  return a.x < b.x;
}

bool compareY(KDCoord a, KDCoord b) {
  return a.y < b.y;
}

KDTree::KDTree(const double *x, const double *y, const int num) {
  n = 0;
  vector<KDCoord> points;
  for(int i = 0; i < num; i++) {
    if(!isnan(x[i]) && !isnan(y[i])) {
      n++;
      KDCoord pt;
      pt.x = x[i];
      pt.y = y[i];
      points.push_back(pt);
    }
  }

  construct_tree(points.begin(), points.end(), 0);
}

int KDTree::construct_tree(vector<KDCoord>::iterator pts_start, vector<KDCoord>::iterator pts_end, int axis) {
  if(axis == 0) {
    sort(pts_start, pts_end, compareX);
  } else {
    sort(pts_start, pts_end, compareY);
  }

  // Create node with median point for this axis
  vector<KDCoord>::iterator median = pts_start + (pts_end - pts_start) / 2;
  KDNode node;
  node.l = -1;
  node.r = -1;
  node.coord = *median;
  nodes.push_back(node);
  int node_ind = nodes.size() - 1;
  int left_child  = -1;
  int right_child = -1;

  // Recursive calls
  if(pts_end - pts_start > 1) {
    if(median - pts_start >= 1)
      left_child = construct_tree(pts_start, median, (axis + 1) % 2);
    if(pts_end - (median + 1) >= 1)
      right_child = construct_tree(median + 1, pts_end, (axis + 1) % 2);
  }

  // Set children after recursive calls (to prevent seg fault caused by the vector being reallocated)
  nodes[node_ind].l = left_child;
  nodes[node_ind].r = right_child;

  return node_ind;
}
