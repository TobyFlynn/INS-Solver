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

double KDTree::closest_point(const double x, const double y) {
  double closest_distance = -1.0;
  int closest_ind = -1;
  int current_ind = 0;
  int axis = 0;

  nearest_neighbour(x, y, current_ind, axis, closest_ind, closest_distance);
  return closest_distance;
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

void KDTree::nearest_neighbour(const double x, const double y, int current_ind, int axis, int &closest_ind, double &closest_distance) {
  if(current_ind == -1) return;
  double distance = (x - nodes[current_ind].coord.x) * (x - nodes[current_ind].coord.x)
                     + (y - nodes[current_ind].coord.y) * (y - nodes[current_ind].coord.y);
  // Check if leaf node
  if(nodes[current_ind].l == -1 && nodes[current_ind].r == -1) {
    // Check if first leaf node
    if(closest_ind == -1) {
      closest_ind = current_ind;
      closest_distance = distance;
      return;
    }

    // Check if this leaf node is closer than current best
    if(distance < closest_distance) {
      closest_ind = current_ind;
      closest_distance = distance;
    }
    return;
  }

  // If not leaf node, search tree as if inserting on the way down
  if(axis == 0) {
    if(x < nodes[current_ind].coord.x) {
      nearest_neighbour(x, y, nodes[current_ind].l, (axis + 1) % 2, closest_ind, closest_distance);
    } else {
      nearest_neighbour(x, y, nodes[current_ind].r, (axis + 1) % 2, closest_ind, closest_distance);
    }
  } else {
    if(y < nodes[current_ind].coord.y) {
      nearest_neighbour(x, y, nodes[current_ind].l, (axis + 1) % 2, closest_ind, closest_distance);
    } else {
      nearest_neighbour(x, y, nodes[current_ind].r, (axis + 1) % 2, closest_ind, closest_distance);
    }
  }

  // Check if current node is better than current best
  if(distance < closest_distance || closest_ind == -1) {
    closest_ind = current_ind;
    closest_distance = distance;
  }

  // Check if there could be a closer point on the other side of the plane
  if(axis == 0) {
    double plane_distance = (x - nodes[current_ind].coord.x) * (x - nodes[current_ind].coord.x);
    if(plane_distance > closest_distance) {
      // No possible better node on other side of this node
      return;
    } else {
      // Potentially better node on other side so search there
      if(x < nodes[current_ind].coord.x) {
        nearest_neighbour(x, y, nodes[current_ind].r, (axis + 1) % 2, closest_ind, closest_distance);
      } else {
        nearest_neighbour(x, y, nodes[current_ind].l, (axis + 1) % 2, closest_ind, closest_distance);
      }
    }
  } else {
    double plane_distance = (y - nodes[current_ind].coord.y) * (y - nodes[current_ind].coord.y);
    if(plane_distance > closest_distance) {
      // No possible better node on other side of this node
      return;
    } else {
      // Potentially better node on other side so search there
      if(y < nodes[current_ind].coord.y) {
        nearest_neighbour(x, y, nodes[current_ind].r, (axis + 1) % 2, closest_ind, closest_distance);
      } else {
        nearest_neighbour(x, y, nodes[current_ind].l, (axis + 1) % 2, closest_ind, closest_distance);
      }
    }
  }
}
