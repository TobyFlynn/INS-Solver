#include "kd_tree.h"

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
  return a.x < b.x;
}

bool compareY(KDCoord a, KDCoord b) {
  return a.y < b.y;
}

KDTree::KDTree(const double *x, const double *y, const int num, 
               DGMesh *mesh, op_dat s) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  n = 0;
  vector<KDCoord> points;
  for(int i = 0; i < num; i++) {
    if(!isnan(x[i]) && !isnan(y[i])) {
      n++;
      KDCoord pt;
      pt.x = x[i];
      pt.y = y[i];
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
  int blockLenghts[]  = {2, 2};
  MPI_Aint displacements[] = {offsetof(KDCoord, x), offsetof(KDCoord, poly)};
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
  construct_tree(points.begin(), points.end(), 0);
}

KDCoord KDTree::closest_point(const double x, const double y) {
  double closest_distance = -1.0;
  int closest_ind = -1;
  int current_ind = 0;
  int axis = 0;
  KDCoord result;

  nearest_neighbour(x, y, current_ind, axis, closest_ind, closest_distance, result);

  return result;
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

void KDTree::nearest_neighbour(const double x, const double y, int current_ind, int axis, int &closest_ind, double &closest_distance, KDCoord &pt) {
  if(current_ind == -1) return;
  double distance = (x - nodes[current_ind].coord.x) * (x - nodes[current_ind].coord.x)
                     + (y - nodes[current_ind].coord.y) * (y - nodes[current_ind].coord.y);
  // Check if leaf node
  if(nodes[current_ind].l == -1 && nodes[current_ind].r == -1) {
    // Check if first leaf node
    if(closest_ind == -1) {
      closest_ind = current_ind;
      closest_distance = distance;
      pt = nodes[current_ind].coord;
      return;
    }

    // Check if this leaf node is closer than current best
    if(distance < closest_distance) {
      closest_ind = current_ind;
      closest_distance = distance;
      pt = nodes[current_ind].coord;
    }
    return;
  }

  // If not leaf node, search tree as if inserting on the way down
  if(axis == 0) {
    if(x < nodes[current_ind].coord.x) {
      nearest_neighbour(x, y, nodes[current_ind].l, (axis + 1) % 2, closest_ind, closest_distance, pt);
    } else {
      nearest_neighbour(x, y, nodes[current_ind].r, (axis + 1) % 2, closest_ind, closest_distance, pt);
    }
  } else {
    if(y < nodes[current_ind].coord.y) {
      nearest_neighbour(x, y, nodes[current_ind].l, (axis + 1) % 2, closest_ind, closest_distance, pt);
    } else {
      nearest_neighbour(x, y, nodes[current_ind].r, (axis + 1) % 2, closest_ind, closest_distance, pt);
    }
  }

  // Check if current node is better than current best
  if(distance < closest_distance || closest_ind == -1) {
    closest_ind = current_ind;
    closest_distance = distance;
    pt = nodes[current_ind].coord;
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
        nearest_neighbour(x, y, nodes[current_ind].r, (axis + 1) % 2, closest_ind, closest_distance, pt);
      } else {
        nearest_neighbour(x, y, nodes[current_ind].l, (axis + 1) % 2, closest_ind, closest_distance, pt);
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
        nearest_neighbour(x, y, nodes[current_ind].r, (axis + 1) % 2, closest_ind, closest_distance, pt);
      } else {
        nearest_neighbour(x, y, nodes[current_ind].l, (axis + 1) % 2, closest_ind, closest_distance, pt);
      }
    }
  }
}

std::set<int> KDTree::cell_inds(vector<KDCoord> &points) {
  std::set<int> result;
  for(int i = 0; i < points.size(); i++) {
    result.insert(points[i].poly);
  }
  return result;
}

void KDTree::construct_polys(vector<KDCoord> &points, DGMesh *mesh, op_dat s) {
  timer->startTimer("LS - Construct Poly Approx");
  // Get cell inds that require polynomial approximations
  std::set<int> cellInds = cell_inds(points);

  map<int,set<int>> stencils = PolyApprox::get_stencils(cellInds, mesh->edge2cells);

  const double *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const double *y_ptr = getOP2PtrHost(mesh->y, OP_READ);
  const double *s_ptr = getOP2PtrHost(s, OP_READ);

  // Populate map
  int i = 0;
  for(auto it = cellInds.begin(); it != cellInds.end(); it++) {
    set<int> stencil = stencils.at(*it);
    PolyApprox p(*it, stencil, x_ptr, y_ptr, s_ptr);
    polys.push_back(p);
    cell2polyMap.insert({*it, i});
    i++;
  }

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHost(s, OP_READ, s_ptr);

  timer->endTimer("LS - Construct Poly Approx");
}

vector<PolyApprox> KDTree::get_polys() {
  return polys;
}

void KDTree::update_poly_inds(vector<KDCoord> &points) {
  for(auto &point : points) {
    point.poly = cell2polyMap.at(point.poly);
  }
}
