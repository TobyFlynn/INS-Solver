#include "kd_tree_mpi.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstddef>

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

KDTreeMPI::KDTreeMPI(const double *x, const double *y, const int num, 
                     DGMesh *mesh, op_dat s) {
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
      points.push_back(pt);
    }
  }

  // Construct cell to poly map for all these points
  construct_polys(points, mesh, s);
  update_poly_inds(points);

  // Construct local complete tree
  construct_tree(points.begin(), points.end(), false, 0);
}

vector<vector<KDCoord>::iterator> KDTreeMPI::local_search(const int num_pts, const double *x, const double *y) {
  vector<vector<KDCoord>::iterator> local_closest;
  for(int i = 0; i < num_pts; i++) {
    double closest_distance = std::numeric_limits<double>::max();
    vector<KDCoord>::iterator res = points.end();
    int current_ind = 0;
    nearest_neighbour(x[i], y[i], current_ind, res, closest_distance);
    local_closest.push_back(res);
  }
  return local_closest;
}

vector<vector<KDCoord>::iterator> KDTreeMPI::local_search(const int num_pts, const double *pts) {
  vector<vector<KDCoord>::iterator> local_closest;
  for(int i = 0; i < num_pts; i++) {
    double closest_distance = std::numeric_limits<double>::max();
    vector<KDCoord>::iterator res = points.end();
    int current_ind = 0;
    nearest_neighbour(pts[i * 2], pts[i * 2 + 1], current_ind, res, closest_distance);
    local_closest.push_back(res);
  }
  return local_closest;
}

void KDTreeMPI::get_global_bounding_boxes(MPI_Comm *mpi_comm, MPIBB *mpi_bb) {
  double bounding_box[4];
  if(nodes.size() > 0) {
    bounding_box[0] = nodes[0].x_min;
    bounding_box[1] = nodes[0].x_max;
    bounding_box[2] = nodes[0].y_min;
    bounding_box[3] = nodes[0].y_max;
  } else {
    bounding_box[0] = NAN;
    bounding_box[1] = NAN;
    bounding_box[2] = NAN;
    bounding_box[3] = NAN;
  }

  int comm_size;
  MPI_Comm_size(*mpi_comm, &comm_size);

  double *bb_recv_buf = (double *)calloc(4 * comm_size, sizeof(double));
  MPI_Allgather(bounding_box, 4, MPI_DOUBLE, bb_recv_buf, 4, MPI_DOUBLE, *mpi_comm);

  for(int i = 0; i < comm_size; i++) {
    mpi_bb[i].x_min = bb_recv_buf[i * 4];
    mpi_bb[i].x_max = bb_recv_buf[i * 4 + 1];
    mpi_bb[i].y_min = bb_recv_buf[i * 4 + 2];
    mpi_bb[i].y_max = bb_recv_buf[i * 4 + 3];
  }
  free(bb_recv_buf);
}

void KDTreeMPI::round1_get_pts_to_send_to_ranks(const int num_pts, const double *x, const double *y,
                                                const MPIBB *mpi_bb, const int Reinit_comm_size,
                                                int *num_pts_to_send, map<int,vector<int>> &rank_to_pt_inds) {
  // Work out which points to send to which processor
  for(int i = 0; i < num_pts; i++) {
    // For each point check which bounding box is closest
    double closest_dist = std::numeric_limits<double>::max();
    int closest_rank = -1;
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(!isnan(mpi_bb[rank].x_min)) {
        double dist = bb_sqr_dist(mpi_bb[rank], x[i], y[i]);
        if(dist < closest_dist) {
          closest_dist = dist;
          closest_rank = rank;
        }
      }
    }
    // Add 1 to num of points to send to closest rank
    num_pts_to_send[closest_rank]++;
    // Add this pt index to the list of pts to send to closest rank
    rank_to_pt_inds.at(closest_rank).push_back(i);
  }
}

void KDTreeMPI::round1_prepare_send_recv(const int num_pts, const double *x, const double *y, 
                                         const int Reinit_comm_size, map<int,vector<int>> &rank_to_pt_inds,
                                         int *num_pts_to_recv, int *send_inds, int *recv_inds, 
                                         double **pts_to_send, double **pts_to_recv, vector<int> &pt_send_rcv_map,
                                         int &num_remote_pts) {
  if(nodes.size() == 0) {
    *pts_to_send = (double *)calloc(num_pts * 2, sizeof(double));
    double *pts_send = *pts_to_send;
    int ind = 0;
    // Iterate over each rank
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      send_inds[rank] = ind * 2;
      // Iterate over all points to send to this specific rank
      for(const auto &pt : rank_to_pt_inds.at(rank)) {
        pts_send[ind * 2]     = x[pt];
        pts_send[ind * 2 + 1] = y[pt];
        pt_send_rcv_map.push_back(pt);
        ind++;
      }
    }
  } else {
    num_remote_pts = 0;
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      recv_inds[rank] = num_remote_pts * 2;
      num_remote_pts += num_pts_to_recv[rank];
    }
    *pts_to_recv = (double *)calloc(num_remote_pts * 2, sizeof(double));
  }
}

void KDTreeMPI::round1_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm, 
                             int *num_pts_to_send, int *num_pts_to_recv, double *pts_to_send, double *pts_to_recv, 
                             int *send_inds, int *recv_inds, MPI_Request *requests) {
  if(nodes.size() == 0) {
    // The sends
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_send[rank] != 0) {
        double *send_buf = &pts_to_send[send_inds[rank]];
        MPI_Isend(send_buf, 2 * num_pts_to_send[rank], MPI_DOUBLE, rank, 0, *mpi_comm, &requests[rank]);
      }
    }
  } else {
    // The receives
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_recv[rank] != 0) {
        double *recv_buf = &pts_to_recv[recv_inds[rank]];
        MPI_Irecv(recv_buf, 2 * num_pts_to_recv[rank], MPI_DOUBLE, rank, 0, *mpi_comm, &requests[rank]);
      }
    }
  }
}

void KDTreeMPI::round1_wait_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Request *requests, 
                                  int *num_pts_to_send, int *num_pts_to_recv) {
  if(nodes.size() == 0) {
    // The sends
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_send[rank] != 0) {
        MPI_Wait(&requests[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&requests[rank]);
      }
    }
  } else {
    // The receives
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_recv[rank] != 0) {
        MPI_Wait(&requests[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&requests[rank]);
      }
    }
  }
}

void KDTreeMPI::round1_send_results(const int Reinit_comm_rank, const int Reinit_comm_size, const int num_pts, const int num_remote_pts,
                                    int *num_pts_to_send, int *num_pts_to_recv, int *send_inds, int *recv_inds,
                                    vector<vector<KDCoord>::iterator> &remote_closest, MPIKDResponse **response,
                                    MPI_Datatype *mpi_type, MPI_Comm *mpi_comm) {
  MPI_Request requests[Reinit_comm_size];
  if(nodes.size() == 0) {
    *response = (MPIKDResponse *)calloc(num_pts, sizeof(MPIKDResponse));
    MPIKDResponse *resp = *response;
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_send[rank] != 0) {
        MPIKDResponse *recv_buf = &resp[send_inds[rank] / 2];
        MPI_Irecv(recv_buf, num_pts_to_send[rank], *mpi_type, rank, 0, *mpi_comm, &requests[rank]);
      }
    }
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_send[rank] != 0) {
        MPI_Wait(&requests[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&requests[rank]);
      }
    }
  } else {
    *response = (MPIKDResponse *)calloc(num_remote_pts, sizeof(MPIKDResponse));
    MPIKDResponse *resp = *response;
    for(int i = 0; i < remote_closest.size(); i++) {
      resp[i].x = remote_closest[i]->x;
      resp[i].y = remote_closest[i]->y;
      resp[i].poly = remote_closest[i]->poly;
    }
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_recv[rank] != 0) {
        MPIKDResponse *send_buf = &resp[recv_inds[rank] / 2];
        MPI_Isend(send_buf, num_pts_to_recv[rank], *mpi_type, rank, 0, *mpi_comm, &requests[rank]);
      }
    }
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_recv[rank] != 0) {
        MPI_Wait(&requests[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&requests[rank]);
      }
    }
  }
}

void KDTreeMPI::closest_point(const int num_pts, const double *x, const double *y) {
  // 2) Get comm of ranks that contain nodes to be reinitialised (num_pts != 0)
  MPI_Comm Reinit_comm;
  int Reinit_comm_size, Reinit_comm_rank;
  if(num_pts != 0.0) {
    MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &Reinit_comm);
    MPI_Comm_size(Reinit_comm, &Reinit_comm_size);
    MPI_Comm_rank(Reinit_comm, &Reinit_comm_rank);
  } else {
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, 0, &Reinit_comm);
  }

  if(num_pts == 0)
    return;

  // 3) Share bounding box of each MPI rank
  MPIBB mpi_bb[Reinit_comm_size];
  get_global_bounding_boxes(&Reinit_comm, mpi_bb);

  // 4) Non-blocking communication of ranks that have no local k-d tree
  //    saying which process they are sending their points too
  int num_pts_to_send[Reinit_comm_size];
  for(int i = 0; i < Reinit_comm_size; i++) {
    num_pts_to_send[i] = 0;
  }
  map<int,vector<int>> rank_to_pt_inds;
  for(int i = 0; i < Reinit_comm_size; i++) {
    vector<int> empty_vec;
    rank_to_pt_inds.insert({i, empty_vec});
  }
  if(nodes.size() == 0) {
    round1_get_pts_to_send_to_ranks(num_pts, x, y, mpi_bb, Reinit_comm_size,
                                    num_pts_to_send, rank_to_pt_inds);
  }

  // Tell each process how many points to expect to receive from each process
  int num_pts_to_recv[Reinit_comm_size];
  MPI_Alltoall(num_pts_to_send, 1, MPI_INT, num_pts_to_recv, 1, MPI_INT, Reinit_comm);

  // Pack points to send to each process
  double *pts_to_send;
  double *pts_to_recv;
  int send_inds[Reinit_comm_size];
  int recv_inds[Reinit_comm_size];
  vector<int> pt_send_rcv_map;
  int num_remote_pts;
  round1_prepare_send_recv(num_pts, x, y, Reinit_comm_size, rank_to_pt_inds,
                           num_pts_to_recv,send_inds, recv_inds, 
                           &pts_to_send, &pts_to_recv, pt_send_rcv_map, num_remote_pts);

  // The actual non blocking comms
  MPI_Request requests[Reinit_comm_size];
  round1_comms(Reinit_comm_rank, Reinit_comm_size, &Reinit_comm, num_pts_to_send, num_pts_to_recv,
               pts_to_send, pts_to_recv, send_inds, recv_inds, requests);

  // 5) Ranks with a local k-d tree perform local k-d tree searches,
  //    overlaping with step 4
  vector<vector<KDCoord>::iterator> local_closest;
  if(nodes.size() > 0) {
    local_closest = local_search(num_pts, x, y);
  }

  // 6) Wait on step 4 and then process these nodes on the relevant ranks
  round1_wait_comms(Reinit_comm_rank, Reinit_comm_size, requests, 
                    num_pts_to_send, num_pts_to_recv);

  // Now search local k-d tree using these nodes
  vector<vector<KDCoord>::iterator> remote_closest;
  if(nodes.size() > 0) {
    remote_closest = local_search(num_remote_pts, pts_to_recv);
  }

  // 7) Send back the results of step 6
  MPIKDResponse *response;
  int blockLenghts[]  = {2, 1};
  MPI_Aint displacements[] = {offsetof(MPIKDResponse, x), offsetof(MPIKDResponse, poly)};
  MPI_Datatype types[] = {MPI_DOUBLE, MPI_INT};
  MPI_Datatype MPI_MPIKDResponse_Type;
  MPI_Type_create_struct(2, blockLenghts, displacements, types, &MPI_MPIKDResponse_Type);
  MPI_Type_commit(&MPI_MPIKDResponse_Type);
  round1_send_results(Reinit_comm_rank, Reinit_comm_size, num_pts, num_remote_pts,
                      num_pts_to_send, num_pts_to_recv, send_inds, recv_inds,
                      remote_closest, &response, &MPI_MPIKDResponse_Type, &Reinit_comm);

  // 8) Blocking communication of points to ranks that could potentially 
  //    contain closer points

  // Work out which points could have closer points
  vector<QueryPt> queryPoints;
  if(nodes.size() > 0) {
    for(int i = 0; i < num_pts; i++) {
      QueryPt qp;
      qp.ind = i;
      qp.x = x[i];
      qp.y = y[i];
      qp.closest_x = local_closest[i]->x;
      qp.closest_y = local_closest[i]->y;
      qp.poly = local_closest[i]->poly;
      qp.closest_rank = Reinit_comm_rank;
      qp.lockedin = check_if_locked_in(qp, Reinit_comm_size, mpi_bb);
      queryPoints.push_back(qp);
    }
  } else {
    int currentRank = 0;
    int num_pts_to_send_acc[Reinit_comm_size];
    num_pts_to_send_acc[0] = num_pts_to_send[0];
    for(int i = 1; i < Reinit_comm_size; i++) {
      num_pts_to_send_acc[i] = num_pts_to_send_acc[i - 1] + num_pts_to_send[i];
    }
    for(int i = 0; i < num_pts; i++) {
      while(num_pts_to_send_acc[currentRank] <= i) {
        currentRank++;
      }
      QueryPt qp;
      qp.ind = pt_send_rcv_map[i];
      qp.x = x[qp.ind];
      qp.y = y[qp.ind];
      qp.closest_x = response[i].x;
      qp.closest_y = response[i].y;
      qp.poly = response[i].poly;
      qp.closest_rank = currentRank;
      qp.lockedin = check_if_locked_in(qp, Reinit_comm_size, mpi_bb);
      queryPoints.push_back(qp);
    }
  }

  // Get non locked in points that will need to be sent
  vector<QueryPt> nonLockedIn;
  for(auto &qp : queryPoints) {
    if(!qp.lockedin) {
      nonLockedIn.push_back(qp);
    }
  }
  cout << "Number of non locked in query points: " << nonLockedIn.size() << endl;

  // 9) Do search for other ranks' points

  // 10) Return results

  // 11) Combine remote and local searches to get actual closest point
  //     and associated polynomial approximation

  // 12) Free MPI stuff
  MPI_Type_free(&MPI_MPIKDResponse_Type);
  free(response);
  if(nodes.size() == 0) {
    free(pts_to_send);
  } else {
    free(pts_to_recv);
  }
  MPI_Comm_free(&Reinit_comm);
}

int KDTreeMPI::construct_tree(vector<KDCoord>::iterator pts_start, vector<KDCoord>::iterator pts_end, bool has_transformed, int level) {
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

double KDTreeMPI::bb_sqr_dist(const int node_ind, const double x, const double y) {
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

double KDTreeMPI::bb_sqr_dist(const MPIBB bb, const double x, const double y) {
  double sqr_dist = 0.0;
  if(x < bb.x_min)
    sqr_dist += (x - bb.x_min) * (x - bb.x_min);
  else if(x > bb.x_max)
    sqr_dist += (x - bb.x_max) * (x - bb.x_max);
  
  if(y < bb.y_min)
    sqr_dist += (y - bb.y_min) * (y - bb.y_min);
  else if(y > bb.y_max)
    sqr_dist += (y - bb.y_max) * (y - bb.y_max);
  
  return sqr_dist;
}

bool KDTreeMPI::check_if_locked_in(QueryPt &qp, const int num_ranks, const MPIBB *bb) {
  double sqr_dist = (qp.x - qp.closest_x) * (qp.x - qp.closest_x) + (qp.y - qp.closest_y) * (qp.y - qp.closest_y);
  for(int rank = 0; rank < num_ranks; rank++) {
    if(rank != qp.closest_rank && !isnan(bb[rank].x_min)) {
      double bb_dist = bb_sqr_dist(bb[rank], qp.x, qp.y);
      if(bb_dist < sqr_dist) return false;
    }
  }
  return true;
}

void KDTreeMPI::nearest_neighbour(double x, double y, int current_ind, vector<KDCoord>::iterator &closest_pt, double &closest_distance) {
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

std::set<int> KDTreeMPI::cell_inds(vector<KDCoord> &points) {
  std::set<int> result;
  for(int i = 0; i < points.size(); i++) {
    result.insert(points[i].poly);
  }
  return result;
}

void KDTreeMPI::construct_polys(vector<KDCoord> &points, DGMesh *mesh, op_dat s) {
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

vector<PolyApprox> KDTreeMPI::get_polys() {
  return polys;
}

void KDTreeMPI::update_poly_inds(vector<KDCoord> &points) {
  for(auto &point : points) {
    point.poly = cell2polyMap.at(point.poly);
  }
}
