#include "ls_utils/3d/kd_tree.h"

#include "timing.h"
#include "op2_utils.h"

extern Timing *timer;

using namespace std;

bool compareX(KDCoord a, KDCoord b) {
  return a.x_rot < b.x_rot;
}

bool compareY(KDCoord a, KDCoord b) {
  return a.y_rot < b.y_rot;
}

bool compareZ(KDCoord a, KDCoord b) {
  return a.z_rot < b.z_rot;
}

KDTree3D::KDTree3D(DGMesh3D *m) {
  mesh = m;
}

void KDTree3D::reset() {
  points.clear();
  nodes.clear();
  cell2polyMap.clear();
  polys.clear();
  n = 0;
}

void KDTree3D::pre_build_setup(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                               const int num, op_dat s) {
  n = 0;
  for(int i = 0; i < num; i++) {
    if(!isnan(x[i]) && !isnan(y[i]) && !isnan(z[i])) {
      n++;
      KDCoord pt;
      pt.x = x[i];
      pt.y = y[i];
      pt.z = z[i];
      pt.x_rot = x[i];
      pt.y_rot = y[i];
      pt.z_rot = z[i];
      pt.poly = i / LS_SAMPLE_NP;
      points.push_back(pt);
    }
  }

  // Construct cell to poly map for all these points
  construct_polys(points, s);
  update_poly_inds(points);
}

void KDTree3D::build_tree(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                          const int num, op_dat s) {
  pre_build_setup(x, y, z, num, s);

  if(points.size() == 0) {
    empty = true;
    return;
  }
  empty = false;

  timer->startTimer("K-D Tree - Construct Tree");
  construct_tree(points.begin(), points.end(), false, 0);
  timer->endTimer("K-D Tree - Construct Tree");
}

KDCoord KDTree3D::closest_point(DG_FP x, DG_FP y, DG_FP z) {
  DG_FP closest_distance = std::numeric_limits<DG_FP>::max();
  vector<KDCoord>::iterator res = points.end();
  int current_ind = 0;

  nearest_neighbour(x, y, z, current_ind, res, closest_distance);

  return *res;
}

int KDTree3D::construct_tree(vector<KDCoord>::iterator pts_start, vector<KDCoord>::iterator pts_end, bool has_transformed, int level) {
  KDNode node;
  node.l = -1;
  node.r = -1;

  // Bounding box and mean
  node.x_min = pts_start->x_rot;
  node.x_max = pts_start->x_rot;
  node.y_min = pts_start->y_rot;
  node.y_max = pts_start->y_rot;
  node.z_min = pts_start->z_rot;
  node.z_max = pts_start->z_rot;
  DG_FP x_avg = pts_start->x_rot;
  DG_FP y_avg = pts_start->y_rot;
  DG_FP z_avg = pts_start->z_rot;
  for(auto it = pts_start + 1; it != pts_end; it++) {
    x_avg += it->x_rot;
    y_avg += it->y_rot;
    z_avg += it->z_rot;
    if(node.x_min > it->x_rot) node.x_min = it->x_rot;
    if(node.y_min > it->y_rot) node.y_min = it->y_rot;
    if(node.z_min > it->z_rot) node.z_min = it->z_rot;
    if(node.x_max < it->x_rot) node.x_max = it->x_rot;
    if(node.y_max < it->y_rot) node.y_max = it->y_rot;
    if(node.z_max < it->z_rot) node.z_max = it->z_rot;
  }
  x_avg /= (DG_FP)(pts_end - pts_start);
  y_avg /= (DG_FP)(pts_end - pts_start);
  z_avg /= (DG_FP)(pts_end - pts_start);

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
  if(node.x_max - node.x_min < node.y_max - node.y_min
     || node.x_max - node.x_min < node.z_max - node.z_min) {
    if(node.z_max - node.z_min < node.y_max - node.y_min) {
      axis = 1;
    } else {
      axis = 2;
    }
  }

  // Do rotational transform if necessary
  bool transform = !has_transformed && level > 5 && pts_end - pts_start >= leaf_size * 4;
  transform = false;
  if(transform) {
    DG_FP hole_radius_sqr = 0.0;
    if(axis == 0) {
      hole_radius_sqr = 0.05 * (node.x_max - node.x_min) * 0.05 * (node.x_max - node.x_min);
    } else if(axis == 1) {
      hole_radius_sqr = 0.05 * (node.y_max - node.y_min) * 0.05 * (node.y_max - node.y_min);
    } else {
      hole_radius_sqr = 0.05 * (node.z_max - node.z_min) * 0.05 * (node.z_max - node.z_min);
    }

    arma::vec normal(3);
    normal(0) = 1.0;
    normal(1) = 0.0;
    normal(2) = 0.0;

    for(auto it = pts_start; it != pts_end; it++) {
      DG_FP x_tmp = it->x_rot - x_avg;
      DG_FP y_tmp = it->y_rot - y_avg;
      DG_FP z_tmp = it->z_rot - z_avg;
      DG_FP msqr = x_tmp * x_tmp + y_tmp * y_tmp + z_tmp * z_tmp;
      if(msqr > hole_radius_sqr) {
        DG_FP tmp_dot = x_tmp * normal(0) + y_tmp * normal(1) + z_tmp * normal(2);
        normal(0) -= x_tmp * tmp_dot / msqr;
        normal(1) -= y_tmp * tmp_dot / msqr;
        normal(2) -= z_tmp * tmp_dot / msqr;
      }
    }

    DG_FP msqr = normal(0) * normal(0) + normal(1) * normal(1) + normal(2) * normal(2);
    if(msqr == 0.0) {
      normal(0) = 1.0;
    } else {
      normal(0) /= sqrt(msqr);
      normal(1) /= sqrt(msqr);
      normal(2) /= sqrt(msqr);
    }

    DG_FP min_alpha = pts_start->x_rot * normal(0) + pts_start->y_rot * normal(1);
    DG_FP max_alpha = min_alpha;
    for(auto it = pts_start + 1; it != pts_end; it++) {
      DG_FP alpha = it->x_rot * normal[0] + it->y_rot * normal[1] + it->z_rot * normal[2];
      if(alpha > max_alpha) max_alpha = alpha;
      if(alpha < min_alpha) min_alpha = alpha;
    }

    DG_FP max_extent = node.x_max - node.x_min;
    if(axis == 1)
      max_extent = node.y_max - node.y_min;
    if(axis == 2)
      max_extent = node.z_max - node.z_min;

    if(max_alpha - min_alpha < 0.1 * max_extent) {
      // Calculate orthonormal basis via Householder matrix
      arma::mat axes(3, 3);
      int j = 0;
      if(fabs(normal(1)) < fabs(normal(0)) || fabs(normal(2)) < fabs(normal(0))) {
        if(fabs(normal(1)) < fabs(normal(2))) {
          j = 1;
        } else {
          j = 2;
        }
      }
      arma::vec u = normal;
      u(j) -= 1.0;
      DG_FP u_norm = sqrt(u(0) * u(0) + u(1) * u(1) + u(2) * u(2));
      u = u / u_norm;
      for(int dim = 0; dim < 3; dim++) {
        for(int i = 0; i < 3; i++) {
          axes(dim, i) = (dim == i ? 1.0 : 0.0) - 2.0 * u(dim) * u(i);
        }
      }

      // Apply coord transformation
      DG_FP alpha_x = axes(0,0) * pts_start->x_rot + axes(0,1) * pts_start->y_rot + axes(0,2) * pts_start->z_rot;
      DG_FP alpha_y = axes(1,0) * pts_start->x_rot + axes(1,1) * pts_start->y_rot + axes(1,2) * pts_start->z_rot;
      DG_FP alpha_z = axes(2,0) * pts_start->x_rot + axes(2,1) * pts_start->y_rot + axes(2,2) * pts_start->z_rot;
      pts_start->x_rot = alpha_x;
      pts_start->y_rot = alpha_y;
      pts_start->z_rot = alpha_z;
      DG_FP b_min_x = alpha_x;
      DG_FP b_max_x = alpha_x;
      DG_FP b_min_y = alpha_y;
      DG_FP b_max_y = alpha_y;
      DG_FP b_min_z = alpha_z;
      DG_FP b_max_z = alpha_z;

      for(auto it = pts_start + 1; it != pts_end; it++) {
        alpha_x = axes(0,0) * it->x_rot + axes(0,1) * it->y_rot + axes(0,2) * it->z_rot;
        alpha_y = axes(1,0) * it->x_rot + axes(1,1) * it->y_rot + axes(1,2) * it->z_rot;
        alpha_z = axes(2,0) * it->x_rot + axes(2,1) * it->y_rot + axes(2,2) * it->z_rot;
        it->x_rot = alpha_x;
        it->y_rot = alpha_y;
        it->z_rot = alpha_y;
        if(alpha_x < b_min_x) b_min_x = alpha_x;
        if(alpha_x > b_max_x) b_max_x = alpha_x;
        if(alpha_y < b_min_y) b_min_y = alpha_y;
        if(alpha_y > b_max_y) b_max_y = alpha_y;
        if(alpha_z < b_min_z) b_min_z = alpha_z;
        if(alpha_z > b_max_z) b_max_z = alpha_z;
      }

      node.rot = axes;
      axis = 0;
      if(b_max_x - b_min_x < b_max_y - b_min_y || b_max_x - b_min_x < b_max_z - b_min_z) {
        if(b_max_z - b_min_z < b_max_y - b_min_y) {
          axis = 1;
        } else {
          axis = 2;
        }
      }
      has_transformed = true;
    }
  }

  if(axis == 0) {
    sort(pts_start, pts_end, compareX);
  } else if(axis == 1) {
    sort(pts_start, pts_end, compareY);
  } else {
    sort(pts_start, pts_end, compareZ);
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
    if(pts_end - (median + 1) >= 1)
      right_child = construct_tree(median + 1, pts_end, has_transformed, level + 1);
  }

  // Set children after recursive calls (to prevent seg fault caused by the vector being reallocated)
  nodes[node_ind].l = left_child;
  nodes[node_ind].r = right_child;

  return node_ind;
}

DG_FP KDTree3D::bb_sqr_dist(const int node_ind, const DG_FP x, const DG_FP y, const DG_FP z) {
  DG_FP sqr_dist = 0.0;
  if(x < nodes[node_ind].x_min)
    sqr_dist += (x - nodes[node_ind].x_min) * (x - nodes[node_ind].x_min);
  else if(x > nodes[node_ind].x_max)
    sqr_dist += (x - nodes[node_ind].x_max) * (x - nodes[node_ind].x_max);

  if(y < nodes[node_ind].y_min)
    sqr_dist += (y - nodes[node_ind].y_min) * (y - nodes[node_ind].y_min);
  else if(y > nodes[node_ind].y_max)
    sqr_dist += (y - nodes[node_ind].y_max) * (y - nodes[node_ind].y_max);

  if(z < nodes[node_ind].z_min)
    sqr_dist += (z - nodes[node_ind].z_min) * (z - nodes[node_ind].z_min);
  else if(z > nodes[node_ind].z_max)
    sqr_dist += (z - nodes[node_ind].z_max) * (z - nodes[node_ind].z_max);

  return sqr_dist;
}

void KDTree3D::nearest_neighbour(DG_FP x, DG_FP y, DG_FP z, int current_ind, vector<KDCoord>::iterator &closest_pt, DG_FP &closest_distance) {
  if(nodes[current_ind].l == -1 && nodes[current_ind].r == -1) {
    // Leaf node
    for(auto it = nodes[current_ind].start; it != nodes[current_ind].end; it++) {
      DG_FP sqr_dist = (it->x_rot - x) * (it->x_rot - x) + (it->y_rot - y) * (it->y_rot - y) + (it->z_rot - z) * (it->z_rot - z);
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
    DG_FP new_x = nodes[current_ind].rot(0,0) * x + nodes[current_ind].rot(0,1) * y + nodes[current_ind].rot(0,2) * z;
    DG_FP new_y = nodes[current_ind].rot(1,0) * x + nodes[current_ind].rot(1,1) * y + nodes[current_ind].rot(1,2) * z;
    DG_FP new_z = nodes[current_ind].rot(2,0) * x + nodes[current_ind].rot(2,1) * y + nodes[current_ind].rot(2,2) * z;
    x = new_x;
    y = new_y;
    z = new_z;
  }

  DG_FP dist_l = bb_sqr_dist(nodes[current_ind].l, x, y, z);
  DG_FP dist_r = bb_sqr_dist(nodes[current_ind].r, x, y, z);

  if(dist_l < dist_r) {
    if(dist_l < closest_distance) {
      nearest_neighbour(x, y, z, nodes[current_ind].l, closest_pt, closest_distance);
      if(dist_r < closest_distance) {
        nearest_neighbour(x, y, z, nodes[current_ind].r, closest_pt, closest_distance);
      }
    }
  } else {
    if(dist_r < closest_distance) {
      nearest_neighbour(x, y, z, nodes[current_ind].r, closest_pt, closest_distance);
      if(dist_l < closest_distance) {
        nearest_neighbour(x, y, z, nodes[current_ind].l, closest_pt, closest_distance);
      }
    }
  }
}

std::set<int> KDTree3D::cell_inds(vector<KDCoord> &points) {
  std::set<int> result;
  for(int i = 0; i < points.size(); i++) {
    result.insert(points[i].poly);
  }
  return result;
}

void KDTree3D::construct_polys(vector<KDCoord> &points, op_dat s) {
  timer->startTimer("LS - Construct Poly Approx");
  // Get cell inds that require polynomial approximations
  std::set<int> cellInds = cell_inds(points);

  map<int,set<int>> stencils = PolyApprox3D::get_stencils(cellInds, mesh->face2cells);

  const DG_FP *x_ptr = getOP2PtrHostHE(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHostHE(mesh->y, OP_READ);
  const DG_FP *z_ptr = getOP2PtrHostHE(mesh->z, OP_READ);
  const DG_FP *s_ptr = getOP2PtrHostHE(s, OP_READ);

  // Populate map
  int i = 0;
  for(auto it = cellInds.begin(); it != cellInds.end(); it++) {
    set<int> stencil = stencils.at(*it);
    PolyApprox3D p(*it, stencil, x_ptr, y_ptr, z_ptr, s_ptr);
    polys.push_back(p);
    cell2polyMap.insert({*it, i});
    i++;
  }

  releaseOP2PtrHostHE(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHostHE(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHostHE(mesh->z, OP_READ, z_ptr);
  releaseOP2PtrHostHE(s, OP_READ, s_ptr);

  timer->endTimer("LS - Construct Poly Approx");
}

vector<PolyApprox3D> KDTree3D::get_polys() {
  return polys;
}

void KDTree3D::update_poly_inds(vector<KDCoord> &points) {
  for(auto &point : points) {
    point.poly = cell2polyMap.at(point.poly);
  }
}
