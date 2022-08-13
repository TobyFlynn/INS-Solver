#include "kd_tree_mpi_naive.h"

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

KDTreeMPINaive::KDTreeMPINaive(const double *x, const double *y, const int num, 
               DGMesh *mesh, op_dat s) {
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

  // Create datatype for KDCoord
  int blockLenghts[]  = {4, 2};
  MPI_Aint displacements[] = {offsetof(KDCoord, x_rot), offsetof(KDCoord, poly)};
  MPI_Datatype types[] = {MPI_DOUBLE, MPI_INT};
  MPI_Datatype MPI_KDCoord_Type;

  MPI_Type_create_struct(2, blockLenghts, displacements, types, &MPI_KDCoord_Type);
  MPI_Type_commit(&MPI_KDCoord_Type);

  // Get info necessary for all gather v of KDCoords
  int *sizes = (int *)malloc(comm_size * sizeof(int));
  int *disp  = (int *)malloc(comm_size * sizeof(int));
  MPI_Allgather(&n, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);
  int curDisp = 0;
  for(int i = 0; i < comm_size; i++) {
    disp[i] = curDisp;
    curDisp += sizes[i];
  }
  int numCoords_g = curDisp;

  // All gather v of KDCoords
  KDCoord *recvbuf = (KDCoord *)malloc(numCoords_g * sizeof(KDCoord));
  MPI_Allgatherv(points.data(), n, MPI_KDCoord_Type, recvbuf, sizes, disp, MPI_KDCoord_Type, MPI_COMM_WORLD);

  // Pack PolyApprox data
  int coeffPerPoly = PolyApprox::num_coeff();
  int polyStride = 2 + coeffPerPoly;
  int numPoly = polys.size();
  double *polySendBuf = (double *)malloc(polyStride * numPoly * sizeof(double));
  int polyInd = 0;
  for(auto &poly : polys) {
    int ind = polyInd * polyStride;
    double x, y;
    poly.get_offsets(x, y);
    polySendBuf[ind]     = x;
    polySendBuf[ind + 1] = y;

    for(int i = 0 ; i < coeffPerPoly; i++) {
      polySendBuf[ind + i + 2] = poly.get_coeff(i);
    }

    polyInd++;
  }

  // Get data needed for all gather v of PolyApprox
  MPI_Allgather(&numPoly, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);
  curDisp = 0;
  for(int i = 0; i < comm_size; i++) {
    disp[i] = curDisp;
    curDisp += sizes[i];
  }
  int numPoly_g = curDisp;

  // All gather v of PolyApprox
  MPI_Datatype MPI_Poly_type;
  MPI_Type_contiguous(polyStride, MPI_DOUBLE, &MPI_Poly_type);
  MPI_Type_commit(&MPI_Poly_type);
  
  double *recvbufPoly = (double *)malloc(numPoly_g * polyStride * sizeof(double));
  MPI_Allgatherv(polySendBuf, numPoly, MPI_Poly_type, recvbufPoly, sizes, disp, MPI_Poly_type, MPI_COMM_WORLD);

  // Reconstruct points vector and poly map
  points.clear();
  for(int i = 0; i < numCoords_g; i++) {
    KDCoord currentCoord = recvbuf[i];
    // Update poly ind now that polys from all ranks have been gathered
    currentCoord.poly = disp[currentCoord.rank] + currentCoord.poly;
    points.push_back(currentCoord);
  }

  polys.clear();
  for(int i = 0; i < numPoly_g; i++) {
    int ind = i * polyStride;
    double x = recvbufPoly[ind];
    double y = recvbufPoly[ind + 1];
    std::vector<double> co;
    for(int j = 0; j < coeffPerPoly; j++) {
      co.push_back(recvbufPoly[ind + 2 + j]);
    }
    PolyApprox p(co, x, y);
    polys.push_back(p);
  }

  // Clean up
  free(recvbufPoly);
  MPI_Type_free(&MPI_Poly_type);
  free(polySendBuf);
  free(disp);
  free(sizes);
  MPI_Type_free(&MPI_KDCoord_Type);
  free(recvbuf);

  // Construct local complete tree
  construct_tree(points.begin(), points.end(), false, 0);
}

KDCoord KDTreeMPINaive::closest_point(const double x, const double y) {
  double closest_distance = std::numeric_limits<double>::max();
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
  double x_avg = pts_start->x_rot;
  double y_avg = pts_start->y_rot;
  for(auto it = pts_start + 1; it != pts_end; it++) {
    x_avg += it->x_rot;
    y_avg += it->y_rot;
    if(node.x_min > it->x_rot) node.x_min = it->x_rot;
    if(node.y_min > it->y_rot) node.y_min = it->y_rot;
    if(node.x_max < it->x_rot) node.x_max = it->x_rot;
    if(node.y_max < it->y_rot) node.y_max = it->y_rot;
  }
  x_avg /= (double)(pts_end - pts_start);
  y_avg /= (double)(pts_end - pts_start);

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
    double hole_radius_sqr = 0.0;
    if(axis == 0) {
      hole_radius_sqr = 0.05 * (node.x_max - node.x_min) * 0.05 * (node.x_max - node.x_min);
    } else {
      hole_radius_sqr = 0.05 * (node.y_max - node.y_min) * 0.05 * (node.y_max - node.y_min);
    }

    arma::vec normal(2);
    normal(0) = 1.0;
    normal(1) = 0.0;

    for(auto it = pts_start; it != pts_end; it++) {
      double x_tmp = it->x_rot - x_avg;
      double y_tmp = it->y_rot - y_avg;
      double msqr = x_tmp * x_tmp + y_tmp * y_tmp;
      if(msqr > hole_radius_sqr) {
        double tmp_dot = x_tmp * normal(0) + y_tmp * normal(1);
        normal(0) -= x_tmp * tmp_dot / msqr;
        normal(1) -= y_tmp * tmp_dot / msqr;
      }
    }

    double msqr = normal(0) * normal(0) + normal(1) * normal(1);
    if(msqr == 0.0) {
      normal(0) = 1.0;
    } else {
      normal(0) /= sqrt(msqr);
      normal(1) /= sqrt(msqr);
    }

    double min_alpha = pts_start->x_rot * normal(0) + pts_start->y_rot * normal(1);
    double max_alpha = min_alpha;
    for(auto it = pts_start + 1; it != pts_end; it++) {
      double alpha = it->x_rot * normal[0] + it->y_rot * normal[1];
      if(alpha > max_alpha) max_alpha = alpha;
      if(alpha < min_alpha) min_alpha = alpha;
    }

    double max_extent = node.x_max - node.x_min;
    if(axis == 1)
      max_extent = node.y_max - node.y_min;
    if(max_alpha - min_alpha < 0.1 * max_extent) {
      // Calculate orthonormal basis via Householder matrix
      arma::mat axes(2, 2);
      int j = fabs(normal(0)) < fabs(normal(1)) ? 0 : 1;
      arma::vec u = normal;
      u(j) -= 1.0;
      double u_norm = sqrt(u(0) * u(0) + u(1) * u(1));
      u = u / u_norm;
      for(int dim = 0; dim < 2; dim++) {
        for(int i = 0; i < 2; i++) {
          axes(dim, i) = (dim == i ? 1.0 : 0.0) - 2.0 * u(dim) * u(i);
        }
      }

      // Apply coord transformation
      double alpha_x = axes(0,0) * pts_start->x_rot + axes(0,1) * pts_start->y_rot;
      double alpha_y = axes(1,0) * pts_start->x_rot + axes(1,1) * pts_start->y_rot;
      pts_start->x_rot = alpha_x;
      pts_start->y_rot = alpha_y;
      double b_min_x = alpha_x;
      double b_max_x = alpha_x;
      double b_min_y = alpha_y;
      double b_max_y = alpha_y;

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

double KDTreeMPINaive::bb_sqr_dist(const int node_ind, const double x, const double y) {
  double sqr_dist = 0.0;
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

void KDTreeMPINaive::nearest_neighbour(double x, double y, int current_ind, vector<KDCoord>::iterator &closest_pt, double &closest_distance) {
  if(nodes[current_ind].l == -1 && nodes[current_ind].r == -1) {
    // Leaf node
    for(auto it = nodes[current_ind].start; it != nodes[current_ind].end; it++) {
      double sqr_dist = (it->x_rot - x) * (it->x_rot - x) + (it->y_rot - y) * (it->y_rot - y);
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
    double new_x = nodes[current_ind].rot(0,0) * x + nodes[current_ind].rot(0,1) * y;
    double new_y = nodes[current_ind].rot(1,0) * x + nodes[current_ind].rot(1,1) * y;
    x = new_x;
    y = new_y;
  }

  double dist_l = bb_sqr_dist(nodes[current_ind].l, x, y);
  double dist_r = bb_sqr_dist(nodes[current_ind].r, x, y);

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

void KDTreeMPINaive::construct_polys(vector<KDCoord> &points, DGMesh *mesh, op_dat s) {
  timer->startTimer("LS - Construct Poly Approx");
  // Get cell inds that require polynomial approximations
  std::set<int> cellInds = cell_inds(points);

  map<int,set<int>> stencils = PolyApprox::get_stencils(cellInds, mesh->edge2cells);

  const double *x_ptr = getOP2PtrHostMap(mesh->x, mesh->edge2cells, OP_READ);
  const double *y_ptr = getOP2PtrHostMap(mesh->y, mesh->edge2cells, OP_READ);
  const double *s_ptr = getOP2PtrHostMap(s, mesh->edge2cells, OP_READ);

  // Populate map
  int i = 0;
  for(auto it = cellInds.begin(); it != cellInds.end(); it++) {
    set<int> stencil = stencils.at(*it);
    PolyApprox p(*it, stencil, x_ptr, y_ptr, s_ptr);
    polys.push_back(p);
    cell2polyMap.insert({*it, i});
    i++;
  }

  releaseOP2PtrHostMap(mesh->x, mesh->edge2cells, OP_READ, x_ptr);
  releaseOP2PtrHostMap(mesh->y, mesh->edge2cells, OP_READ, y_ptr);
  releaseOP2PtrHostMap(s, mesh->edge2cells, OP_READ, s_ptr);

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
