#include "ls_utils/2d/kd_tree_mpi_naive.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstddef>

#include "mpi.h"

#include "timing.h"
#include "utils.h"

extern Timing *timer;

using namespace std;

bool compareX(KDCoord a, KDCoord b) {
  return a.x_rot < b.x_rot;
}

bool compareY(KDCoord a, KDCoord b) {
  return a.y_rot < b.y_rot;
}

KDTreeMPINaive::KDTreeMPINaive(const DG_FP *x, const DG_FP *y, const int num,
               DGMesh2D *mesh, op_dat s) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  n = 0;
  for(int i = 0; i < num; i++) {
    if(!isnan(x[i]) && !isnan(y[i])) {
      n++;
      KDCoord pt;
      pt.x = x[i];
      pt.y = y[i];
      pt.x_rot = x[i];
      pt.y_rot = y[i];
      pt.poly = i / LS_SAMPLE_NP;
      pt.rank = rank;
      points.push_back(pt);
    }
  }

  // Construct cell to poly map for all these points
  construct_polys(points, mesh, s);
  update_poly_inds(points);

  // Initial MPI MVP is to just gather all and then construct full local kd-tree
  // Will obviously improve after getting this working as this involves a
  // ridiculous amount of comms
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  int num_pts_per_rank[comm_size];
  MPI_Allgather(&n, 1, MPI_INT, num_pts_per_rank, 1, MPI_INT, MPI_COMM_WORLD);
  int num_pts_g = 0;
  int num_pts_per_rank_disp[comm_size];
  for(int i = 0; i < comm_size; i++) {
    num_pts_per_rank_disp[i] = num_pts_g;
    num_pts_g += num_pts_per_rank[i];
  }

  int num_polys_per_rank[comm_size];
  int num_polys = polys.size();
  MPI_Allgather(&num_polys, 1, MPI_INT, num_polys_per_rank, 1, MPI_INT, MPI_COMM_WORLD);
  int num_polys_g = 0;
  int num_polys_coeff_per_rank[comm_size];
  int num_polys_coeff_per_rank_disp[comm_size];
  for(int i = 0; i < comm_size; i++) {
    num_polys_coeff_per_rank_disp[i] = num_polys_g * PolyApprox::num_coeff();
    num_polys_g += num_polys_per_rank[i];
    num_polys_coeff_per_rank[i] = num_polys_per_rank[i] * PolyApprox::num_coeff();
  }

  vector<DG_FP> pts_x_snd(n);
  vector<DG_FP> pts_y_snd(n);
  vector<int> pts_poly_snd(n);
  for(int i = 0; i < n; i++) {
    pts_x_snd[i] = points[i].x;
    pts_y_snd[i] = points[i].y;
    pts_poly_snd[i] = points[i].poly;
  }

  vector<DG_FP> poly_coeff_snd(polys.size() * PolyApprox::num_coeff());
  for(int i = 0; i < polys.size(); i++) {
    int ind = i * PolyApprox::num_coeff();
    for(int j = 0; j < PolyApprox::num_coeff(); j++) {
      poly_coeff_snd[ind + j] = polys[i].get_coeff(j);
    }
  }

  vector<DG_FP> pts_x_rcv(num_pts_g);
  vector<DG_FP> pts_y_rcv(num_pts_g);
  vector<int> pts_poly_rcv(num_pts_g);
  vector<DG_FP> poly_coeff_rcv(num_polys_g * PolyApprox::num_coeff());

  MPI_Allgatherv(pts_x_snd.data(), n, DG_MPI_FP, pts_x_rcv.data(), num_pts_per_rank, num_pts_per_rank_disp, DG_MPI_FP, MPI_COMM_WORLD);
  MPI_Allgatherv(pts_y_snd.data(), n, DG_MPI_FP, pts_y_rcv.data(), num_pts_per_rank, num_pts_per_rank_disp, DG_MPI_FP, MPI_COMM_WORLD);
  MPI_Allgatherv(pts_poly_snd.data(), n, MPI_INT, pts_poly_rcv.data(), num_pts_per_rank, num_pts_per_rank_disp, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgatherv(poly_coeff_snd.data(), num_polys * PolyApprox::num_coeff(), DG_MPI_FP, poly_coeff_rcv.data(), num_polys_coeff_per_rank, num_polys_coeff_per_rank_disp, DG_MPI_FP, MPI_COMM_WORLD);

  vector<KDCoord> newPts;
  int curr_rank = 0;
  for(int i = 0; i < num_pts_g; i++) {
    while(curr_rank < comm_size - 1 && i >= num_pts_per_rank_disp[curr_rank + 1])
      curr_rank++;
    KDCoord kc;
    kc.x = pts_x_rcv[i];
    kc.y = pts_y_rcv[i];
    kc.x_rot = pts_x_rcv[i];
    kc.y_rot = pts_y_rcv[i];
    kc.poly = (num_polys_coeff_per_rank_disp[curr_rank] / PolyApprox::num_coeff()) + pts_poly_rcv[i];
    newPts.push_back(kc);
  }

  vector<PolyApprox> newPolys;
  for(int i = 0; i < num_polys_g; i++) {
    vector<DG_FP> c;
    int ind = i * PolyApprox::num_coeff();
    for(int j = 0; j < PolyApprox::num_coeff(); j++) {
      c.push_back(poly_coeff_rcv[ind + j]);
    }
    PolyApprox p(c, 0.0, 0.0);
    newPolys.push_back(p);
  }

  points = newPts;
  polys = newPolys;

  // Construct local complete tree
  construct_tree(points.begin(), points.end(), false, 0);
}

KDCoord KDTreeMPINaive::closest_point(const DG_FP x, const DG_FP y) {
  DG_FP closest_distance = std::numeric_limits<DG_FP>::max();
  vector<KDCoord>::iterator res = points.end();
  int current_ind = 0;

  nearest_neighbour(x, y, current_ind, res, closest_distance);

  return *res;
}

int KDTreeMPINaive::construct_tree(vector<KDCoord>::iterator pts_start, vector<KDCoord>::iterator pts_end, bool has_transformed, int level) {
  KDNode node;
  node.l = -1;
  node.r = -1;

  // Bounding box and mean
  node.x_min = pts_start->x_rot;
  node.x_max = pts_start->x_rot;
  node.y_min = pts_start->y_rot;
  node.y_max = pts_start->y_rot;
  DG_FP x_avg = pts_start->x_rot;
  DG_FP y_avg = pts_start->y_rot;
  for(auto it = pts_start + 1; it != pts_end; it++) {
    x_avg += it->x_rot;
    y_avg += it->y_rot;
    if(node.x_min > it->x_rot) node.x_min = it->x_rot;
    if(node.y_min > it->y_rot) node.y_min = it->y_rot;
    if(node.x_max < it->x_rot) node.x_max = it->x_rot;
    if(node.y_max < it->y_rot) node.y_max = it->y_rot;
  }
  x_avg /= (DG_FP)(pts_end - pts_start);
  y_avg /= (DG_FP)(pts_end - pts_start);

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
  if(node.x_max - node.x_min < node.y_max - node.y_min)
    axis = 1;

  // Do rotational transform if necessary
  bool transform = !has_transformed && level > 5 && pts_end - pts_start >= leaf_size * 4;
  if(transform) {
    DG_FP hole_radius_sqr = 0.0;
    if(axis == 0) {
      hole_radius_sqr = 0.05 * (node.x_max - node.x_min) * 0.05 * (node.x_max - node.x_min);
    } else {
      hole_radius_sqr = 0.05 * (node.y_max - node.y_min) * 0.05 * (node.y_max - node.y_min);
    }

    arma::vec normal(2);
    normal(0) = 1.0;
    normal(1) = 0.0;

    for(auto it = pts_start; it != pts_end; it++) {
      DG_FP x_tmp = it->x_rot - x_avg;
      DG_FP y_tmp = it->y_rot - y_avg;
      DG_FP msqr = x_tmp * x_tmp + y_tmp * y_tmp;
      if(msqr > hole_radius_sqr) {
        DG_FP tmp_dot = x_tmp * normal(0) + y_tmp * normal(1);
        normal(0) -= x_tmp * tmp_dot / msqr;
        normal(1) -= y_tmp * tmp_dot / msqr;
      }
    }

    DG_FP msqr = normal(0) * normal(0) + normal(1) * normal(1);
    if(msqr == 0.0) {
      normal(0) = 1.0;
    } else {
      normal(0) /= sqrt(msqr);
      normal(1) /= sqrt(msqr);
    }

    DG_FP min_alpha = pts_start->x_rot * normal(0) + pts_start->y_rot * normal(1);
    DG_FP max_alpha = min_alpha;
    for(auto it = pts_start + 1; it != pts_end; it++) {
      DG_FP alpha = it->x_rot * normal[0] + it->y_rot * normal[1];
      if(alpha > max_alpha) max_alpha = alpha;
      if(alpha < min_alpha) min_alpha = alpha;
    }

    DG_FP max_extent = node.x_max - node.x_min;
    if(axis == 1)
      max_extent = node.y_max - node.y_min;
    if(max_alpha - min_alpha < 0.1 * max_extent) {
      // Calculate orthonormal basis via Householder matrix
      arma::mat axes(2, 2);
      int j = fabs(normal(0)) < fabs(normal(1)) ? 0 : 1;
      arma::vec u = normal;
      u(j) -= 1.0;
      DG_FP u_norm = sqrt(u(0) * u(0) + u(1) * u(1));
      u = u / u_norm;
      for(int dim = 0; dim < 2; dim++) {
        for(int i = 0; i < 2; i++) {
          axes(dim, i) = (dim == i ? 1.0 : 0.0) - 2.0 * u(dim) * u(i);
        }
      }

      // Apply coord transformation
      DG_FP alpha_x = axes(0,0) * pts_start->x_rot + axes(0,1) * pts_start->y_rot;
      DG_FP alpha_y = axes(1,0) * pts_start->x_rot + axes(1,1) * pts_start->y_rot;
      pts_start->x_rot = alpha_x;
      pts_start->y_rot = alpha_y;
      DG_FP b_min_x = alpha_x;
      DG_FP b_max_x = alpha_x;
      DG_FP b_min_y = alpha_y;
      DG_FP b_max_y = alpha_y;

      for(auto it = pts_start + 1; it != pts_end; it++) {
        alpha_x = axes(0,0) * it->x_rot + axes(0,1) * it->y_rot;
        alpha_y = axes(1,0) * it->x_rot + axes(1,1) * it->y_rot;
        it->x_rot = alpha_x;
        it->y_rot = alpha_y;
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
    if(pts_end - (median + 1) >= 1)
      right_child = construct_tree(median + 1, pts_end, has_transformed, level + 1);
  }

  // Set children after recursive calls (to prevent seg fault caused by the vector being reallocated)
  nodes[node_ind].l = left_child;
  nodes[node_ind].r = right_child;

  return node_ind;
}

DG_FP KDTreeMPINaive::bb_sqr_dist(const int node_ind, const DG_FP x, const DG_FP y) {
  DG_FP sqr_dist = 0.0;
  if(x < nodes[node_ind].x_min)
    sqr_dist += (x - nodes[node_ind].x_min) * (x - nodes[node_ind].x_min);
  else if(x > nodes[node_ind].x_max)
    sqr_dist += (x - nodes[node_ind].x_max) * (x - nodes[node_ind].x_max);

  if(y < nodes[node_ind].y_min)
    sqr_dist += (y - nodes[node_ind].y_min) * (y - nodes[node_ind].y_min);
  else if(y > nodes[node_ind].y_max)
    sqr_dist += (y - nodes[node_ind].y_max) * (y - nodes[node_ind].y_max);

  return sqr_dist;
}

void KDTreeMPINaive::nearest_neighbour(DG_FP x, DG_FP y, int current_ind, vector<KDCoord>::iterator &closest_pt, DG_FP &closest_distance) {
  if(nodes[current_ind].l == -1 && nodes[current_ind].r == -1) {
    // Leaf node
    for(auto it = nodes[current_ind].start; it != nodes[current_ind].end; it++) {
      DG_FP sqr_dist = (it->x_rot - x) * (it->x_rot - x) + (it->y_rot - y) * (it->y_rot - y);
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
    DG_FP new_x = nodes[current_ind].rot(0,0) * x + nodes[current_ind].rot(0,1) * y;
    DG_FP new_y = nodes[current_ind].rot(1,0) * x + nodes[current_ind].rot(1,1) * y;
    x = new_x;
    y = new_y;
  }

  DG_FP dist_l = bb_sqr_dist(nodes[current_ind].l, x, y);
  DG_FP dist_r = bb_sqr_dist(nodes[current_ind].r, x, y);

  if(dist_l < dist_r) {
    if(dist_l < closest_distance) {
      nearest_neighbour(x, y, nodes[current_ind].l, closest_pt, closest_distance);
      if(dist_r < closest_distance) {
        nearest_neighbour(x, y, nodes[current_ind].r, closest_pt, closest_distance);
      }
    }
  } else {
    if(dist_r < closest_distance) {
      nearest_neighbour(x, y, nodes[current_ind].r, closest_pt, closest_distance);
      if(dist_l < closest_distance) {
        nearest_neighbour(x, y, nodes[current_ind].l, closest_pt, closest_distance);
      }
    }
  }
}

std::set<int> KDTreeMPINaive::cell_inds(vector<KDCoord> &points) {
  std::set<int> result;
  for(int i = 0; i < points.size(); i++) {
    result.insert(points[i].poly);
  }
  return result;
}

void KDTreeMPINaive::construct_polys(vector<KDCoord> &points, DGMesh2D *mesh, op_dat s) {
  timer->startTimer("LS - Construct Poly Approx");
  // Get cell inds that require polynomial approximations
  std::set<int> cellInds = cell_inds(points);

  map<int,set<int>> stencils = PolyApprox::get_stencils(cellInds, mesh->face2cells);

  const DG_FP *x_ptr = getOP2PtrHostHE(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHostHE(mesh->y, OP_READ);
  const DG_FP *s_ptr = getOP2PtrHostHE(s, OP_READ);

  // Populate map
  int i = 0;
  for(auto it = cellInds.begin(); it != cellInds.end(); it++) {
    set<int> stencil = stencils.at(*it);
    PolyApprox p(*it, stencil, x_ptr, y_ptr, s_ptr);
    polys.push_back(p);
    cell2polyMap.insert({*it, i});
    i++;
  }

  releaseOP2PtrHostHE(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHostHE(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHostHE(s, OP_READ, s_ptr);

  timer->endTimer("LS - Construct Poly Approx");
}

vector<PolyApprox> KDTreeMPINaive::get_polys() {
  return polys;
}

void KDTreeMPINaive::update_poly_inds(vector<KDCoord> &points) {
  for(auto &point : points) {
    point.poly = cell2polyMap.at(point.poly);
  }
}
