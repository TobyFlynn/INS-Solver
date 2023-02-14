#ifndef __INS_KD_TREE_H
#define __INS_KD_TREE_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include <vector>
#include <set>
#include <map>

#include "mpi.h"

#include "dg_mesh/dg_mesh_2d.h"
#include "ls_reinit_poly.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

struct KDCoord {
  DG_FP x_rot;
  DG_FP y_rot;
  DG_FP x;
  DG_FP y;
  int poly;
  int rank;
};

struct KDNode {
  int l;
  int r;
  std::vector<KDCoord>::iterator start;
  std::vector<KDCoord>::iterator end;
  arma::mat rot;
  DG_FP x_min, x_max;
  DG_FP y_min, y_max;
  int axis;
};

struct MPIBB {
  DG_FP x_min;
  DG_FP x_max;
  DG_FP y_min;
  DG_FP y_max;
};

struct MPIKDResponse {
  DG_FP x;
  DG_FP y;
  int poly;
};

struct QueryPt {
  int ind;
  DG_FP x;
  DG_FP y;
  DG_FP closest_x;
  DG_FP closest_y;
  int closest_rank;
  int poly;
  bool lockedin;
  std::vector<int> potential_ranks;
};

class KDTreeMPI {
public:
  KDTreeMPI(const DG_FP *x, const DG_FP *y, const int num, DGMesh2D *mesh, op_dat s);

  void closest_point(const int num_pts, const DG_FP *x, const DG_FP *y, DG_FP *closest_x, DG_FP *closest_y, int *poly_ind);
  std::vector<PolyApprox> get_polys();

private:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, bool has_transformed, int level);
  DG_FP bb_sqr_dist(const int node_ind, const DG_FP x, const DG_FP y);
  DG_FP bb_sqr_dist(const MPIBB bb, const DG_FP x, const DG_FP y);
  bool check_if_locked_in(QueryPt &qp, const int num_ranks, const MPIBB *bb);
  void nearest_neighbour(DG_FP x, DG_FP y, int current_ind, std::vector<KDCoord>::iterator &closest_pt, DG_FP &closest_distance);

  std::set<int> cell_inds(std::vector<KDCoord> &points);
  void construct_polys(std::vector<KDCoord> &points, DGMesh2D *mesh, op_dat s);
  void update_poly_inds(std::vector<KDCoord> &points);

  std::vector<std::vector<KDCoord>::iterator> local_search(const int num_pts, const DG_FP *x, const DG_FP *y);
  std::vector<std::vector<KDCoord>::iterator> local_search(const int num_pts, const DG_FP *pts);

  void get_global_bounding_boxes(MPI_Comm *mpi_comm, MPIBB *mpi_bb);
  void round1_get_pts_to_send_to_ranks(const int num_pts, const DG_FP *x, const DG_FP *y,
                                       const MPIBB *mpi_bb, const int Reinit_comm_size,
                                       int *num_pts_to_send, std::map<int,std::vector<int>> &rank_to_pt_inds);
  void round1_prepare_send_recv(const int num_pts, const DG_FP *x, const DG_FP *y,
                                const int Reinit_comm_size, std::map<int,std::vector<int>> &rank_to_pt_inds,
                                int *num_pts_to_recv, int *send_inds, int *recv_inds,
                                DG_FP **pts_to_send, DG_FP **pts_to_recv, std::vector<int> &pt_send_rcv_map,
                                int &num_remote_pts);
  void round1_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm,
                    int *num_pts_to_send, int *num_pts_to_recv, DG_FP *pts_to_send, DG_FP *pts_to_recv,
                    int *send_inds, int *recv_inds, MPI_Request *requests);
  void round1_wait_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Request *requests,
                         int *num_pts_to_send, int *num_pts_to_recv);
  void round1_send_results(const int Reinit_comm_rank, const int Reinit_comm_size, const int num_pts, const int num_remote_pts,
                           int *num_pts_to_send, int *num_pts_to_recv, int *send_inds, int *recv_inds,
                           std::vector<std::vector<KDCoord>::iterator> &remote_closest, MPIKDResponse **response,
                           MPI_Datatype *mpi_type, MPI_Comm *mpi_comm);
  std::vector<QueryPt> populate_query_pts(const int num_pts, const DG_FP *x, const DG_FP *y, const int Reinit_comm_rank,
                                   const int Reinit_comm_size, MPIBB *mpi_bb, int *num_pts_to_send, MPIKDResponse *response,
                                   std::vector<int> &pt_send_rcv_map, std::vector<std::vector<KDCoord>::iterator> &local_closest);
  void round2_pack_query_pts(const int Reinit_comm_size, int *num_pts_to_send, int *send_inds,
                             std::vector<QueryPt*> &nonLockedIn, DG_FP **round2_pts_to_send,
                             std::vector<QueryPt*> &qp_ptrs);
  void round2_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm, int *num_pts_to_send, int *num_pts_to_recv,
                             int *send_inds, int *recv_inds, DG_FP *round2_pts_to_send, DG_FP **round2_pts_to_recv);
  void round2_results_comm(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm, MPI_Datatype *mpi_type,
                           int *num_pts_to_send, int *num_pts_to_recv, int *send_inds, int *recv_inds, std::vector<std::vector<KDCoord>::iterator> &remote_closest,
                           MPIKDResponse **round2_send_response, MPIKDResponse **round2_recv_response);
  void round2_update_qp(const int Reinit_comm_size, int *num_pts_to_send, std::vector<QueryPt*> &qp_ptrs, MPIKDResponse *round2_recv_response);
  void get_list_of_polys_wanted(const int Reinit_comm_rank, const int Reinit_comm_size, std::vector<QueryPt> &queryPoints,
                                int *poly_recv_inds, int *num_polys_req, std::vector<int> &polys_wanted);
  void send_list_of_poly_inds(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm,
                              int *num_polys_snd, int *num_polys_req, int *poly_send_inds, int *poly_recv_inds,
                              int **poly_list_to_send, std::vector<int> &polys_wanted);
  void send_polys(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm,
                  int *num_polys_snd, int *num_polys_req, int *poly_send_inds, int *poly_recv_inds,
                  int *poly_list_to_send, DG_FP **requested_poly_coeff);
  void update_local_polys(const int Reinit_comm_rank, const int Reinit_comm_size, int *num_polys_req, int *poly_recv_inds,
                          DG_FP *requested_poly_coeff, std::vector<int> &polys_wanted, std::vector<QueryPt> &queryPoints);

  std::vector<KDNode> nodes;
  std::vector<KDCoord> points;
  int n;
  std::map<int,int> cell2polyMap;
  std::vector<PolyApprox> polys;
  static const int leaf_size = 10;
};

#endif
