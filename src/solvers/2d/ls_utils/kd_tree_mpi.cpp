#include "ls_utils/2d/kd_tree_mpi.h"

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

KDTreeMPI::KDTreeMPI(const DG_FP *x, const DG_FP *y, const int num,
                     DGMesh2D *mesh, op_dat s) {
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
  timer->startTimer("K-D Tree - Construct Tree");
  if(n > 0)
    construct_tree(points.begin(), points.end(), false, 0);
  timer->endTimer("K-D Tree - Construct Tree");
}

vector<vector<KDCoord>::iterator> KDTreeMPI::local_search(const int num_pts, const DG_FP *x, const DG_FP *y) {
  vector<vector<KDCoord>::iterator> local_closest;
  for(int i = 0; i < num_pts; i++) {
    DG_FP closest_distance = std::numeric_limits<DG_FP>::max();
    vector<KDCoord>::iterator res = points.end();
    int current_ind = 0;
    nearest_neighbour(x[i], y[i], current_ind, res, closest_distance);
    local_closest.push_back(res);
  }
  return local_closest;
}

vector<vector<KDCoord>::iterator> KDTreeMPI::local_search(const int num_pts, const DG_FP *pts) {
  vector<vector<KDCoord>::iterator> local_closest;
  for(int i = 0; i < num_pts; i++) {
    DG_FP closest_distance = std::numeric_limits<DG_FP>::max();
    vector<KDCoord>::iterator res = points.end();
    int current_ind = 0;
    nearest_neighbour(pts[i * 2], pts[i * 2 + 1], current_ind, res, closest_distance);
    local_closest.push_back(res);
  }
  return local_closest;
}

void KDTreeMPI::get_global_bounding_boxes(MPI_Comm *mpi_comm, MPIBB *mpi_bb) {
  DG_FP bounding_box[4];
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

  DG_FP *bb_recv_buf = (DG_FP *)calloc(4 * comm_size, sizeof(DG_FP));
  MPI_Allgather(bounding_box, 4, DG_MPI_FP, bb_recv_buf, 4, DG_MPI_FP, *mpi_comm);

  for(int i = 0; i < comm_size; i++) {
    mpi_bb[i].x_min = bb_recv_buf[i * 4];
    mpi_bb[i].x_max = bb_recv_buf[i * 4 + 1];
    mpi_bb[i].y_min = bb_recv_buf[i * 4 + 2];
    mpi_bb[i].y_max = bb_recv_buf[i * 4 + 3];
  }
  free(bb_recv_buf);
}

void KDTreeMPI::round1_get_pts_to_send_to_ranks(const int num_pts, const DG_FP *x, const DG_FP *y,
                                                const MPIBB *mpi_bb, const int Reinit_comm_size,
                                                int *num_pts_to_send, map<int,vector<int>> &rank_to_pt_inds) {
  // Work out which points to send to which processor
  for(int i = 0; i < num_pts; i++) {
    // For each point check which bounding box is closest
    DG_FP closest_dist = std::numeric_limits<DG_FP>::max();
    int closest_rank = -1;
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(!isnan(mpi_bb[rank].x_min)) {
        DG_FP dist = bb_sqr_dist(mpi_bb[rank], x[i], y[i]);
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

void KDTreeMPI::round1_prepare_send_recv(const int num_pts, const DG_FP *x, const DG_FP *y,
                                         const int Reinit_comm_size, map<int,vector<int>> &rank_to_pt_inds,
                                         int *num_pts_to_recv, int *send_inds, int *recv_inds,
                                         DG_FP **pts_to_send, DG_FP **pts_to_recv, vector<int> &pt_send_rcv_map,
                                         int &num_remote_pts) {
  if(nodes.size() == 0) {
    *pts_to_send = (DG_FP *)calloc(num_pts * 2, sizeof(DG_FP));
    DG_FP *pts_send = *pts_to_send;
    int ind = 0;
    // Iterate over each rank
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      send_inds[rank] = ind;
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
      recv_inds[rank] = num_remote_pts;
      num_remote_pts += num_pts_to_recv[rank];
    }
    *pts_to_recv = (DG_FP *)calloc(num_remote_pts * 2, sizeof(DG_FP));
  }
}

void KDTreeMPI::round1_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm,
                             int *num_pts_to_send, int *num_pts_to_recv, DG_FP *pts_to_send, DG_FP *pts_to_recv,
                             int *send_inds, int *recv_inds, MPI_Request *requests) {
  if(nodes.size() == 0) {
    // The sends
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_send[rank] != 0) {
        DG_FP *send_buf = &pts_to_send[send_inds[rank] * 2];
        MPI_Isend(send_buf, 2 * num_pts_to_send[rank], DG_MPI_FP, rank, 0, *mpi_comm, &requests[rank]);
      }
    }
  } else {
    // The receives
    for(int rank = 0; rank < Reinit_comm_size; rank++) {
      if(rank != Reinit_comm_rank && num_pts_to_recv[rank] != 0) {
        DG_FP *recv_buf = &pts_to_recv[recv_inds[rank] * 2];
        MPI_Irecv(recv_buf, 2 * num_pts_to_recv[rank], DG_MPI_FP, rank, 0, *mpi_comm, &requests[rank]);
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
        MPIKDResponse *recv_buf = &resp[send_inds[rank]];
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
        MPIKDResponse *send_buf = &resp[recv_inds[rank]];
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

vector<QueryPt> KDTreeMPI::populate_query_pts(const int num_pts, const DG_FP *x, const DG_FP *y, const int Reinit_comm_rank,
                                   const int Reinit_comm_size, MPIBB *mpi_bb, int *num_pts_to_send, MPIKDResponse *response,
                                   vector<int> &pt_send_rcv_map, vector<vector<KDCoord>::iterator> &local_closest) {
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
      while(currentRank < Reinit_comm_size - 1 && i >= num_pts_to_send_acc[currentRank]) {
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
  return queryPoints;
}

void KDTreeMPI::round2_pack_query_pts(const int Reinit_comm_size, int *num_pts_to_send, int *send_inds,
                                      vector<QueryPt*> &nonLockedIn, DG_FP **round2_pts_to_send,
                                      vector<QueryPt*> &qp_ptrs) {
  map<int,vector<int>> rank_to_qp;
  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    vector<int> empty_vec;
    rank_to_qp.insert({rank, empty_vec});
    num_pts_to_send[rank] = 0;
  }
  int total = 0;
  for(int i = 0; i < nonLockedIn.size(); i++) {
    for(auto rank : nonLockedIn[i]->potential_ranks) {
      rank_to_qp.at(rank).push_back(i);
      num_pts_to_send[rank]++;
      total++;
    }
  }
  *round2_pts_to_send = (DG_FP *)calloc(total * 2, sizeof(DG_FP));
  DG_FP *send_ptr = *round2_pts_to_send;
  int count = 0;
  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    send_inds[rank] = count;
    for(int i = 0; i < rank_to_qp.at(rank).size(); i++) {
      int qp_ind = rank_to_qp.at(rank)[i];
      send_ptr[2 * count] = nonLockedIn[qp_ind]->x;
      send_ptr[2 * count + 1] = nonLockedIn[qp_ind]->y;
      qp_ptrs.push_back(nonLockedIn[qp_ind]);
      count++;
    }
  }
}

void KDTreeMPI::round2_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm, int *num_pts_to_send, int *num_pts_to_recv,
                             int *send_inds, int *recv_inds, DG_FP *round2_pts_to_send, DG_FP **round2_pts_to_recv) {
  recv_inds[0] = 0;
  for(int rank = 1; rank < Reinit_comm_size; rank++) {
    recv_inds[rank] = recv_inds[rank - 1] + num_pts_to_recv[rank - 1];
  }
  int total_recv = recv_inds[Reinit_comm_size - 1] + num_pts_to_recv[Reinit_comm_size - 1];
  *round2_pts_to_recv = (DG_FP *)calloc(total_recv * 2, sizeof(DG_FP));
  DG_FP *recv_ptr = *round2_pts_to_recv;

  MPI_Request send_rq[Reinit_comm_size];
  MPI_Request recv_rq[Reinit_comm_size];
  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    if(rank != Reinit_comm_rank) {
      // Send
      if(num_pts_to_send[rank] != 0) {
        DG_FP *send_buf = &round2_pts_to_send[send_inds[rank] * 2];
        MPI_Isend(send_buf, 2 * num_pts_to_send[rank], DG_MPI_FP, rank, 0, *mpi_comm, &send_rq[rank]);
      }
      // Recv
      if(num_pts_to_recv[rank] != 0) {
        DG_FP *recv_buf = &recv_ptr[recv_inds[rank] * 2];
        MPI_Irecv(recv_buf, 2 * num_pts_to_recv[rank], DG_MPI_FP, rank, 0, *mpi_comm, &recv_rq[rank]);
      }
    }
  }

  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    if(rank != Reinit_comm_rank) {
      // Send
      if(num_pts_to_send[rank] != 0) {
        MPI_Wait(&send_rq[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&send_rq[rank]);
      }
      // Recv
      if(num_pts_to_recv[rank] != 0) {
        MPI_Wait(&recv_rq[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&recv_rq[rank]);
      }
    }
  }
}

void KDTreeMPI::round2_results_comm(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm, MPI_Datatype *mpi_type,
                                    int *num_pts_to_send, int *num_pts_to_recv, int *send_inds, int *recv_inds, vector<vector<KDCoord>::iterator> &remote_closest,
                                    MPIKDResponse **round2_send_response, MPIKDResponse **round2_recv_response) {
  *round2_recv_response = (MPIKDResponse *)calloc(send_inds[Reinit_comm_size - 1] + num_pts_to_send[Reinit_comm_size - 1], sizeof(MPIKDResponse));
  *round2_send_response = (MPIKDResponse *)calloc(remote_closest.size(), sizeof(MPIKDResponse));
  MPIKDResponse* send_ptr = *round2_send_response;
  MPIKDResponse* recv_ptr = *round2_recv_response;

  for(int i = 0; i < remote_closest.size(); i++) {
    send_ptr[i].x = remote_closest[i]->x;
    send_ptr[i].y = remote_closest[i]->y;
    send_ptr[i].poly = remote_closest[i]->poly;
  }
  MPI_Request send_rq[Reinit_comm_size];
  MPI_Request recv_rq[Reinit_comm_size];
  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    if(rank != Reinit_comm_rank) {
      // Send
      if(num_pts_to_recv[rank] != 0) {
        MPIKDResponse *send_buf = &send_ptr[recv_inds[rank]];
        MPI_Isend(send_buf, num_pts_to_recv[rank], *mpi_type, rank, 0, *mpi_comm, &send_rq[rank]);
      }
      // Recv
      if(num_pts_to_send[rank] != 0) {
        MPIKDResponse *recv_buf = &recv_ptr[send_inds[rank]];
        MPI_Irecv(recv_buf, num_pts_to_send[rank], *mpi_type, rank, 0, *mpi_comm, &recv_rq[rank]);
      }
    }
  }

  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    if(rank != Reinit_comm_rank) {
      // Send
      if(num_pts_to_recv[rank] != 0) {
        MPI_Wait(&send_rq[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&send_rq[rank]);
      }
      // Recv
      if(num_pts_to_send[rank] != 0) {
        MPI_Wait(&recv_rq[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&recv_rq[rank]);
      }
    }
  }
}

void KDTreeMPI::round2_update_qp(const int Reinit_comm_size, int *num_pts_to_send, vector<QueryPt*> &qp_ptrs, MPIKDResponse *round2_recv_response) {
  int currentRank = 0;
  int num_pts_to_send_acc[Reinit_comm_size];
  num_pts_to_send_acc[0] = num_pts_to_send[0];
  for(int i = 1; i < Reinit_comm_size; i++) {
    num_pts_to_send_acc[i] = num_pts_to_send_acc[i - 1] + num_pts_to_send[i];
  }
  for(int i = 0; i < qp_ptrs.size(); i++) {
    while(currentRank < Reinit_comm_size - 1 && i >= num_pts_to_send_acc[currentRank]) {
      currentRank++;
    }
    DG_FP x = qp_ptrs[i]->x;
    DG_FP y = qp_ptrs[i]->y;
    DG_FP c_x = qp_ptrs[i]->closest_x;
    DG_FP c_y = qp_ptrs[i]->closest_y;
    DG_FP current_dist = (x - c_x) * (x - c_x) + (y - c_y) * (y - c_y);
    DG_FP nc_x = round2_recv_response[i].x;
    DG_FP nc_y = round2_recv_response[i].y;
    DG_FP new_dist = (x - nc_x) * (x - nc_x) + (y - nc_y) * (y - nc_y);

    if(new_dist < current_dist) {
      qp_ptrs[i]->closest_x = nc_x;
      qp_ptrs[i]->closest_y = nc_y;
      qp_ptrs[i]->closest_rank = currentRank;
      qp_ptrs[i]->poly = round2_recv_response[i].poly;
    }
  }
}

void KDTreeMPI::get_list_of_polys_wanted(const int Reinit_comm_rank, const int Reinit_comm_size, vector<QueryPt> &queryPoints,
                                         int *poly_recv_inds, int *num_polys_req, vector<int> &polys_wanted) {
  map<int,set<int>> rank_to_polys;
  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    set<int> empty_set;
    rank_to_polys.insert({rank, empty_set});
  }
  for(const auto &qp : queryPoints) {
    if(qp.closest_rank != Reinit_comm_rank) {
      rank_to_polys.at(qp.closest_rank).insert(qp.poly);
    }
  }
  // Get list of polys wanted
  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    poly_recv_inds[rank] = polys_wanted.size();
    num_polys_req[rank] = rank_to_polys.at(rank).size();
    for(const auto &pl : rank_to_polys.at(rank)) {
      polys_wanted.push_back(pl);
    }
  }
}

void KDTreeMPI::send_list_of_poly_inds(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm,
                                       int *num_polys_snd, int *num_polys_req, int *poly_send_inds, int *poly_recv_inds,
                                       int **poly_list_to_send, vector<int> &polys_wanted) {
  poly_send_inds[0] = 0;
  for(int rank = 1; rank < Reinit_comm_size; rank++) {
    poly_send_inds[rank] = poly_send_inds[rank - 1] + num_polys_snd[rank - 1];
  }
  int total_poly_snd = poly_send_inds[Reinit_comm_size - 1] + num_polys_snd[Reinit_comm_size - 1];
  *poly_list_to_send = (int *)calloc(total_poly_snd, sizeof(int));
  int *polys_send = *poly_list_to_send;

  MPI_Request send_rq[Reinit_comm_size];
  MPI_Request recv_rq[Reinit_comm_size];
  int *poly_list_wanted = polys_wanted.data();
  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    if(rank != Reinit_comm_rank) {
      // Send
      if(num_polys_req[rank] != 0) {
        int *send_buf = &poly_list_wanted[poly_recv_inds[rank]];
        MPI_Isend(send_buf, num_polys_req[rank], MPI_INT, rank, 0, *mpi_comm, &send_rq[rank]);
      }
      // Recv
      if(num_polys_snd[rank] != 0) {
        int *recv_buf = &polys_send[poly_send_inds[rank]];
        MPI_Irecv(recv_buf, num_polys_snd[rank], MPI_INT, rank, 0, *mpi_comm, &recv_rq[rank]);
      }
    }
  }

  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    if(rank != Reinit_comm_rank) {
      // Send
      if(num_polys_req[rank] != 0) {
        MPI_Wait(&send_rq[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&send_rq[rank]);
      }
      // Recv
      if(num_polys_snd[rank] != 0) {
        MPI_Wait(&recv_rq[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&recv_rq[rank]);
      }
    }
  }
}

void KDTreeMPI::send_polys(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm,
                           int *num_polys_snd, int *num_polys_req, int *poly_send_inds, int *poly_recv_inds,
                           int *poly_list_to_send, DG_FP **requested_poly_coeff) {
  int total_poly_snd = poly_send_inds[Reinit_comm_size - 1] + num_polys_snd[Reinit_comm_size - 1];
  DG_FP *poly_coeff_snd = (DG_FP *)calloc(total_poly_snd * PolyApprox::num_coeff(), sizeof(DG_FP));
  int total_poly_req = poly_recv_inds[Reinit_comm_size - 1] + num_polys_req[Reinit_comm_size - 1];
  *requested_poly_coeff = (DG_FP *)calloc(total_poly_req * PolyApprox::num_coeff(), sizeof(DG_FP));
  DG_FP *req_coeff_ptr = *requested_poly_coeff;
  for(int i = 0; i < total_poly_snd; i++) {
    int poly_ind = poly_list_to_send[i];
    int coeff_ind = i * PolyApprox::num_coeff();
    for(int c = 0; c < PolyApprox::num_coeff(); c++) {
      poly_coeff_snd[coeff_ind + c] = polys[poly_ind].get_coeff(c);
    }
  }

  MPI_Request send_rq[Reinit_comm_size];
  MPI_Request recv_rq[Reinit_comm_size];
  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    if(rank != Reinit_comm_rank) {
      // Send
      if(num_polys_snd[rank] != 0) {
        DG_FP *send_buf = &poly_coeff_snd[poly_send_inds[rank] * PolyApprox::num_coeff()];
        MPI_Isend(send_buf, PolyApprox::num_coeff() * num_polys_snd[rank], DG_MPI_FP, rank, 0, *mpi_comm, &send_rq[rank]);
      }
      // Recv
      if(num_polys_req[rank] != 0) {
        DG_FP *recv_buf = &req_coeff_ptr[poly_recv_inds[rank] * PolyApprox::num_coeff()];
        MPI_Irecv(recv_buf, PolyApprox::num_coeff() * num_polys_req[rank], DG_MPI_FP, rank, 0, *mpi_comm, &recv_rq[rank]);
      }
    }
  }

  for(int rank = 0; rank < Reinit_comm_size; rank++) {
    if(rank != Reinit_comm_rank) {
      // Send
      if(num_polys_snd[rank] != 0) {
        MPI_Wait(&send_rq[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&send_rq[rank]);
      }
      // Recv
      if(num_polys_req[rank] != 0) {
        MPI_Wait(&recv_rq[rank], MPI_STATUS_IGNORE);
        MPI_Request_free(&recv_rq[rank]);
      }
    }
  }

  free(poly_coeff_snd);
}

void KDTreeMPI::update_local_polys(const int Reinit_comm_rank, const int Reinit_comm_size, int *num_polys_req, int *poly_recv_inds,
                                   DG_FP *requested_poly_coeff, vector<int> &polys_wanted, vector<QueryPt> &queryPoints) {
  map<pair<int,int>,int> rank_poly_to_poly;
  int total_poly_req = poly_recv_inds[Reinit_comm_size - 1] + num_polys_req[Reinit_comm_size - 1];
  int currentRank = 0;
  int num_polys_req_acc[Reinit_comm_size];
  num_polys_req_acc[0] = num_polys_req[0];
  for(int i = 1; i < Reinit_comm_size; i++) {
    num_polys_req_acc[i] = num_polys_req_acc[i - 1] + num_polys_req[i];
  }
  for(int i = 0; i < total_poly_req; i++) {
    while(currentRank < Reinit_comm_size - 1 && i >= num_polys_req_acc[currentRank]) {
      currentRank++;
    }
    int poly_ind = i * PolyApprox::num_coeff();
    vector<DG_FP> coeff;
    for(int c = 0; c < PolyApprox::num_coeff(); c++) {
      coeff.push_back(requested_poly_coeff[poly_ind + c]);
    }
    // TODO enable non-zero offset
    PolyApprox p(coeff, 0.0, 0.0);
    polys.push_back(p);
    rank_poly_to_poly.insert({{currentRank, polys_wanted[i]}, polys.size() - 1});
  }

  for(auto &qp : queryPoints) {
    if(qp.closest_rank != Reinit_comm_rank) {
      qp.poly = rank_poly_to_poly.at({qp.closest_rank, qp.poly});
      qp.closest_rank = Reinit_comm_rank;
    }
  }
}

void KDTreeMPI::closest_point(const int num_pts, const DG_FP *x, const DG_FP *y, DG_FP *closest_x, DG_FP *closest_y, int *poly_ind) {
  // 2) Get comm of ranks that contain nodes to be reinitialised (num_pts != 0)
  timer->startTimer("K-D Tree - setup comms");
  MPI_Comm Reinit_comm;
  int Reinit_comm_size, Reinit_comm_rank;
  // if(num_pts != 0) {
    MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &Reinit_comm);
    MPI_Comm_size(Reinit_comm, &Reinit_comm_size);
    MPI_Comm_rank(Reinit_comm, &Reinit_comm_rank);
  // } else {
  //   MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, 0, &Reinit_comm);
  // }
  timer->endTimer("K-D Tree - setup comms");

  // if(num_pts == 0)
  //   return;

  // 3) Share bounding box of each MPI rank
  timer->startTimer("K-D Tree - bounding boxes");
  MPIBB mpi_bb[Reinit_comm_size];
  get_global_bounding_boxes(&Reinit_comm, mpi_bb);
  timer->endTimer("K-D Tree - bounding boxes");

  // 4) Non-blocking communication of ranks that have no local k-d tree
  //    saying which process they are sending their points too
  timer->startTimer("K-D Tree - round1");
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
  DG_FP *pts_to_send;
  DG_FP *pts_to_recv;
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
  timer->startTimer("K-D Tree - round1 local search");
  vector<vector<KDCoord>::iterator> local_closest;
  if(nodes.size() > 0) {
    local_closest = local_search(num_pts, x, y);
  }
  timer->endTimer("K-D Tree - round1 local search");

  // 6) Wait on step 4 and then process these nodes on the relevant ranks
  round1_wait_comms(Reinit_comm_rank, Reinit_comm_size, requests,
                    num_pts_to_send, num_pts_to_recv);

  // Now search local k-d tree using these nodes
  timer->endTimer("K-D Tree - round1 remote search");
  vector<vector<KDCoord>::iterator> remote_closest;
  if(nodes.size() > 0) {
    remote_closest = local_search(num_remote_pts, pts_to_recv);
  }
  timer->endTimer("K-D Tree - round1 remote search");

  // 7) Send back the results of step 6
  MPIKDResponse *response;
  int blockLenghts[]  = {2, 1};
  MPI_Aint displacements[] = {offsetof(MPIKDResponse, x), offsetof(MPIKDResponse, poly)};
  MPI_Datatype types[] = {DG_MPI_FP, MPI_INT};
  MPI_Datatype MPI_MPIKDResponse_Type;
  MPI_Type_create_struct(2, blockLenghts, displacements, types, &MPI_MPIKDResponse_Type);
  MPI_Type_commit(&MPI_MPIKDResponse_Type);
  round1_send_results(Reinit_comm_rank, Reinit_comm_size, num_pts, num_remote_pts,
                      num_pts_to_send, num_pts_to_recv, send_inds, recv_inds,
                      remote_closest, &response, &MPI_MPIKDResponse_Type, &Reinit_comm);
  timer->endTimer("K-D Tree - round1");
  // 8) Blocking communication of points to ranks that could potentially
  //    contain closer points

  // Work out which points could have closer points
  timer->startTimer("K-D Tree - populate query pts");
  vector<QueryPt> queryPoints = populate_query_pts(num_pts, x, y, Reinit_comm_rank, Reinit_comm_size,
                                                   mpi_bb, num_pts_to_send, response, pt_send_rcv_map,
                                                   local_closest);
  // Get non locked in points that will need to be sent
  vector<QueryPt*> nonLockedIn;
  for(auto &qp : queryPoints) {
    if(!qp.lockedin) {
      nonLockedIn.push_back(&qp);
    }
  }
  timer->endTimer("K-D Tree - populate query pts");

  // Pack data to send to remote ranks (some pts may be sent to multiple ranks)
  timer->startTimer("K-D Tree - round2");
  DG_FP *round2_pts_to_send;
  vector<QueryPt*> qp_ptrs;
  round2_pack_query_pts(Reinit_comm_size, num_pts_to_send, send_inds,
                        nonLockedIn, &round2_pts_to_send, qp_ptrs);

  // Tell each process how many points to expect to receive from each process
  MPI_Alltoall(num_pts_to_send, 1, MPI_INT, num_pts_to_recv, 1, MPI_INT, Reinit_comm);

  // Round 2 comms (blocking)
  DG_FP *round2_pts_to_recv;
  round2_comms(Reinit_comm_rank, Reinit_comm_size, &Reinit_comm, num_pts_to_send, num_pts_to_recv,
               send_inds, recv_inds, round2_pts_to_send, &round2_pts_to_recv);

  // 9) Do search for other ranks' points
  remote_closest.clear();
  if(nodes.size() > 0) {
    remote_closest = local_search(num_pts_to_recv[Reinit_comm_size - 1] + recv_inds[Reinit_comm_size - 1], round2_pts_to_recv);
  }

  // 10) Return results
  MPIKDResponse *round2_send_response;
  MPIKDResponse *round2_recv_response;
  round2_results_comm(Reinit_comm_rank, Reinit_comm_size, &Reinit_comm, &MPI_MPIKDResponse_Type,
                      num_pts_to_send, num_pts_to_recv, send_inds, recv_inds, remote_closest,
                      &round2_send_response, &round2_recv_response);
  timer->endTimer("K-D Tree - round2");

  // 11) Combine remote and local searches to get actual closest point
  //     and associated polynomial approximation
  timer->startTimer("K-D Tree - work out closest pt");
  round2_update_qp(Reinit_comm_size, num_pts_to_send, qp_ptrs, round2_recv_response);
  timer->endTimer("K-D Tree - work out closest pt");

  // All query points now have the closest point, just need to get the polys from other ranks
  // Get list of polys wanted
  timer->startTimer("K-D Tree - snd/rcv polys");
  int num_polys_req[Reinit_comm_size];
  int poly_recv_inds[Reinit_comm_size];
  vector<int> polys_wanted;
  get_list_of_polys_wanted(Reinit_comm_rank, Reinit_comm_size, queryPoints,
                           poly_recv_inds, num_polys_req, polys_wanted);

  // Tell each process how many polys to send to each process
  int num_polys_snd[Reinit_comm_size];
  MPI_Alltoall(num_polys_req, 1, MPI_INT, num_polys_snd, 1, MPI_INT, Reinit_comm);

  // Send list of required polys
  int poly_send_inds[Reinit_comm_size];
  int *poly_list_to_send;
  send_list_of_poly_inds(Reinit_comm_rank, Reinit_comm_size, &Reinit_comm,
                         num_polys_snd, num_polys_req, poly_send_inds, poly_recv_inds,
                         &poly_list_to_send, polys_wanted);

  // Send requested polys
  DG_FP *requested_poly_coeff;
  send_polys(Reinit_comm_rank, Reinit_comm_size, &Reinit_comm,
             num_polys_snd, num_polys_req, poly_send_inds, poly_recv_inds,
             poly_list_to_send, &requested_poly_coeff);
  timer->endTimer("K-D Tree - snd/rcv polys");

  // Add received polys to local list of polys and update poly indices
  timer->startTimer("K-D Tree - update local polys");
  update_local_polys(Reinit_comm_rank, Reinit_comm_size, num_polys_req, poly_recv_inds,
                     requested_poly_coeff, polys_wanted, queryPoints);
  timer->endTimer("K-D Tree - update local polys");

  // Return results
  for(auto &qp : queryPoints) {
    closest_x[qp.ind] = qp.closest_x;
    closest_y[qp.ind] = qp.closest_y;
    poly_ind[qp.ind] = qp.poly;
  }

  // 12) Free MPI stuff
  free(requested_poly_coeff);
  free(poly_list_to_send);
  free(round2_recv_response);
  free(round2_send_response);
  free(round2_pts_to_recv);
  free(round2_pts_to_send);
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

DG_FP KDTreeMPI::bb_sqr_dist(const int node_ind, const DG_FP x, const DG_FP y) {
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

DG_FP KDTreeMPI::bb_sqr_dist(const MPIBB bb, const DG_FP x, const DG_FP y) {
  DG_FP sqr_dist = 0.0;
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
  DG_FP sqr_dist = (qp.x - qp.closest_x) * (qp.x - qp.closest_x) + (qp.y - qp.closest_y) * (qp.y - qp.closest_y);
  qp.potential_ranks.clear();
  for(int rank = 0; rank < num_ranks; rank++) {
    if(rank != qp.closest_rank && !isnan(bb[rank].x_min)) {
      DG_FP bb_dist = bb_sqr_dist(bb[rank], qp.x, qp.y);
      if(bb_dist < sqr_dist) qp.potential_ranks.push_back(rank);
    }
  }
  return qp.potential_ranks.size() == 0;
}

void KDTreeMPI::nearest_neighbour(DG_FP x, DG_FP y, int current_ind, vector<KDCoord>::iterator &closest_pt, DG_FP &closest_distance) {
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

std::set<int> KDTreeMPI::cell_inds(vector<KDCoord> &points) {
  std::set<int> result;
  for(int i = 0; i < points.size(); i++) {
    result.insert(points[i].poly);
  }
  return result;
}

void KDTreeMPI::construct_polys(vector<KDCoord> &points, DGMesh2D *mesh, op_dat s) {
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

vector<PolyApprox> KDTreeMPI::get_polys() {
  return polys;
}

void KDTreeMPI::update_poly_inds(vector<KDCoord> &points) {
  for(auto &point : points) {
    point.poly = cell2polyMap.at(point.poly);
  }
}
