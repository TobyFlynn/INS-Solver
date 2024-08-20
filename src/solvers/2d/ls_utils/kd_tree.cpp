#include "ls_utils/2d/kd_tree.h"

#include <cmath>
#include <algorithm>
#include <iostream>

#include "timing.h"
#include "op2_utils.h"

extern Timing *timer;

using namespace std;

bool compareX(KDCoord a, KDCoord b) {
  return a.coord_rot[0] < b.coord_rot[0];
}

bool compareY(KDCoord a, KDCoord b) {
  return a.coord_rot[1] < b.coord_rot[1];
}

KDTree::KDTree(DGMesh2D *m) {
  mesh = m;
}

void KDTree::reset() {
  points.clear();
  nodes.clear();
  cell2polyMap.clear();
  polys.clear();
  n = 0;
}

void KDTree::pre_build_setup(const DG_FP *x, const DG_FP *y, const int num, op_dat s) {
  n = 0;
  for(int i = 0; i < num; i++) {
    if(!isnan(x[i])) {
      n++;
    }
  }
  points.resize(n);
  int j = 0;
  for(int i = 0; i < num; i++) {
    if(!isnan(x[i])) {
      points[j].coord[0] = x[i];
      points[j].coord[1] = y[i];
      points[j].coord_rot[0] = x[i];
      points[j].coord_rot[1] = y[i];
      points[j].poly = i / LS_SAMPLE_NP;
      j++;
    }
  }

  // Construct cell to poly map for all these points
  update_poly_inds(points);
}

void KDTree::build_tree(const DG_FP *x, const DG_FP *y, const int num, op_dat s) {
  pre_build_setup(x, y, num, s);

  if(points.size() == 0) {
    empty = true;
    return;
  }
  empty = false;

  timer->startTimer("K-D Tree - Construct Tree");
  construct_tree(points.begin(), points.end(), false, 0);
  timer->endTimer("K-D Tree - Construct Tree");
}

KDCoord KDTree::closest_point(DG_FP x, DG_FP y) {
  DG_FP closest_distance = std::numeric_limits<DG_FP>::max();
  vector<KDCoord>::iterator res = points.end();
  int current_ind = 0;
  DGUtils::Vec<2> coord(x, y);

  nearest_neighbour(coord, current_ind, res, closest_distance);

  return *res;
}

int KDTree::construct_tree(vector<KDCoord>::iterator pts_start, vector<KDCoord>::iterator pts_end, bool has_transformed, int level) {
  KDNode node;
  node.l = -1;
  node.r = -1;

  // Bounding box and mean
  node.min = pts_start->coord_rot;
  node.max = pts_start->coord_rot;
  DGUtils::Vec<2> avg = pts_start->coord_rot;
  for(auto it = pts_start + 1; it != pts_end; it++) {
    avg += it->coord_rot;
    if(node.min[0] > it->coord_rot[0]) node.min[0] = it->coord_rot[0];
    if(node.min[1] > it->coord_rot[1]) node.min[1] = it->coord_rot[1];
    if(node.max[0] < it->coord_rot[0]) node.max[0] = it->coord_rot[0];
    if(node.max[1] < it->coord_rot[1]) node.max[1] = it->coord_rot[1];
  }
  avg /= (DG_FP)(pts_end - pts_start);

  if(pts_end - pts_start <= leaf_size) {
    node.start = pts_start;
    node.end = pts_end;
    nodes.push_back(node);
    int node_ind = nodes.size() - 1;
    return node_ind;
  }

  // Construct splitting node

  // Split across axis with greatest extent
  int axis = 0;
  DGUtils::Vec<2> extent = node.max - node.min;
  DG_FP max_extent = extent[0];
  if(extent[0] < extent[1]) {
    axis = 1;
    max_extent = extent[1];
  }

  // Do rotational transform if necessary
  bool transform = !has_transformed && level > 5 && pts_end - pts_start >= leaf_size * 4;
  if(transform) {
    DG_FP hole_radius_sqr = 0.05 * max_extent * 0.05 * max_extent;

    DGUtils::Vec<2> normal(1.0, 0.0);
    for(auto it = pts_start; it != pts_end; it++) {
      DGUtils::Vec<2> tmp = it->coord_rot - avg;
      DG_FP msqr = tmp.sqr_magnitude();
      if(msqr > hole_radius_sqr) {
        DG_FP dot = tmp.dot(normal);
        normal -= tmp * (dot / msqr);
      }
    }

    DG_FP msqr = normal.sqr_magnitude();
    if(msqr == 0.0) {
      normal[0] = 1.0;
    } else {
      normal /= sqrt(msqr);
    }

    DG_FP min_alpha = pts_start->coord_rot.dot(normal);
    DG_FP max_alpha = min_alpha;
    for(auto it = pts_start + 1; it != pts_end; it++) {
      DG_FP alpha = it->coord_rot.dot(normal);
      if(alpha > max_alpha) max_alpha = alpha;
      if(alpha < min_alpha) min_alpha = alpha;
    }

    if(max_alpha - min_alpha < 0.1 * max_extent) {
      // Calculate orthonormal basis via Householder matrix
      arma::mat axes(2, 2);
      int j = fabs(normal[0]) < fabs(normal[1]) ? 0 : 1;
      DGUtils::Vec<2> u = normal;
      u[j] -= 1.0;
      DG_FP u_norm = u.magnitude();
      u /= u_norm;
      for(int dim = 0; dim < 2; dim++) {
        for(int i = 0; i < 2; i++) {
          axes(dim, i) = (dim == i ? 1.0 : 0.0) - 2.0 * u[dim] * u[i];
        }
      }

      // Apply coord transformation
      DG_FP alpha_x = axes(0,0) * pts_start->coord_rot[0] + axes(0,1) * pts_start->coord_rot[1];
      DG_FP alpha_y = axes(1,0) * pts_start->coord_rot[0] + axes(1,1) * pts_start->coord_rot[1];
      pts_start->coord_rot[0] = alpha_x;
      pts_start->coord_rot[1] = alpha_y;
      DG_FP b_min_x = alpha_x;
      DG_FP b_max_x = alpha_x;
      DG_FP b_min_y = alpha_y;
      DG_FP b_max_y = alpha_y;

      for(auto it = pts_start + 1; it != pts_end; it++) {
        alpha_x = axes(0,0) * it->coord_rot[0] + axes(0,1) * it->coord_rot[1];
        alpha_y = axes(1,0) * it->coord_rot[0] + axes(1,1) * it->coord_rot[1];
        it->coord_rot[0] = alpha_x;
        it->coord_rot[1] = alpha_y;
        if(alpha_x < b_min_x) b_min_x = alpha_x;
        if(alpha_x > b_max_x) b_max_x = alpha_x;
        if(alpha_y < b_min_y) b_min_y = alpha_y;
        if(alpha_y > b_max_y) b_max_y = alpha_y;
      }

      node.rot = axes;
      axis = b_max_x - b_min_x >= b_max_y - b_min_y ? 0 : 1;
      has_transformed = true;
    }
  }

  if(axis == 0) {
    sort(pts_start, pts_end, compareX);
  } else {
    sort(pts_start, pts_end, compareY);
  }

  // Create node with median point for this axis
  vector<KDCoord>::iterator median = pts_start + (pts_end - pts_start) / 2;
  node.axis = axis;
  nodes.push_back(node);
  int node_ind = nodes.size() - 1;
  int left_child  = -1;
  int right_child = -1;

  // Recursive calls
  if(pts_end - pts_start > 1) {
    if(median - pts_start >= 1)
      left_child = construct_tree(pts_start, median, has_transformed, level + 1);
    if(pts_end - median >= 1)
      right_child = construct_tree(median, pts_end, has_transformed, level + 1);
  }

  // Set children after recursive calls (to prevent seg fault caused by the vector being reallocated)
  nodes[node_ind].l = left_child;
  nodes[node_ind].r = right_child;

  return node_ind;
}

DG_FP KDTree::bb_sqr_dist(const int node_ind, DGUtils::Vec<2> &coord) {
  DG_FP sqr_dist = 0.0;
  if(coord[0] < nodes[node_ind].min[0])
    sqr_dist += (coord[0] - nodes[node_ind].min[0]) * (coord[0] - nodes[node_ind].min[0]);
  else if(coord[0] > nodes[node_ind].max[0])
    sqr_dist += (coord[0] - nodes[node_ind].max[0]) * (coord[0] - nodes[node_ind].max[0]);

  if(coord[1] < nodes[node_ind].min[1])
    sqr_dist += (coord[1] - nodes[node_ind].min[1]) * (coord[1] - nodes[node_ind].min[1]);
  else if(coord[1] > nodes[node_ind].max[1])
    sqr_dist += (coord[1] - nodes[node_ind].max[1]) * (coord[1] - nodes[node_ind].max[1]);

  return sqr_dist;
}

void KDTree::nearest_neighbour(DGUtils::Vec<2> coord, int current_ind, vector<KDCoord>::iterator &closest_pt, DG_FP &closest_distance) {
  if(nodes[current_ind].l == -1 && nodes[current_ind].r == -1) {
    // Leaf node
    for(auto it = nodes[current_ind].start; it != nodes[current_ind].end; it++) {
      DG_FP sqr_dist = (it->coord_rot - coord).sqr_magnitude();
      if(closest_distance > sqr_dist) {
        closest_distance = sqr_dist;
        closest_pt = it;
      }
    }
    return;
  }

  // Non leaf node

  // Apply transform
  if(nodes[current_ind].rot.n_elem > 1) {
    DG_FP new_x = nodes[current_ind].rot(0,0) * coord[0] + nodes[current_ind].rot(0,1) * coord[1];
    DG_FP new_y = nodes[current_ind].rot(1,0) * coord[0] + nodes[current_ind].rot(1,1) * coord[1];
    coord[0] = new_x;
    coord[1] = new_y;
  }

  DG_FP dist_l = bb_sqr_dist(nodes[current_ind].l, coord);
  DG_FP dist_r = bb_sqr_dist(nodes[current_ind].r, coord);

  if(dist_l < dist_r) {
    if(dist_l < closest_distance) {
      nearest_neighbour(coord, nodes[current_ind].l, closest_pt, closest_distance);
      if(dist_r < closest_distance) {
        nearest_neighbour(coord, nodes[current_ind].r, closest_pt, closest_distance);
      }
    }
  } else {
    if(dist_r < closest_distance) {
      nearest_neighbour(coord, nodes[current_ind].r, closest_pt, closest_distance);
      if(dist_l < closest_distance) {
        nearest_neighbour(coord, nodes[current_ind].l, closest_pt, closest_distance);
      }
    }
  }
}

void KDTree::set_poly_data(std::vector<PolyApprox> &_polys,
                             std::map<int,int> &_cell2polyMap) {
  polys = _polys;
  cell2polyMap = _cell2polyMap;
}

vector<PolyApprox> KDTree::get_polys() {
  return polys;
}

void KDTree::update_poly_inds(vector<KDCoord> &points) {
  for(auto &point : points) {
    point.poly = cell2polyMap.at(point.poly);
  }
}
