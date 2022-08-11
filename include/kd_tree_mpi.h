#ifndef __INS_KD_TREE_H
#define __INS_KD_TREE_H

#include "op_seq.h"

#include <vector>
#include <set>
#include <map>

#include "mpi.h"

#include "dg_mesh.h"
#include "ls_reinit_poly.h"

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

struct KDCoord {
  double x_rot;
  double y_rot;
  double x;
  double y;
  int poly;
  int rank;
};

struct KDNode {
  int l;
  int r;
  std::vector<KDCoord>::iterator start;
  std::vector<KDCoord>::iterator end;
  arma::mat rot;
  double x_min, x_max;
  double y_min, y_max;
  int axis;
};

struct MPIBB {
  double x_min;
  double x_max;
  double y_min;
  double y_max;
};

struct MPIKDResponse {
  double x;
  double y;
  int poly;
};

struct QueryPt {
  int ind;
  double x;
  double y;
  double closest_x;
  double closest_y;
  int closest_rank;
  int poly;
  bool lockedin;
  std::vector<int> potential_ranks;
};

class KDTreeMPI {
public:
  KDTreeMPI(const double *x, const double *y, const int num, DGMesh *mesh, op_dat s);

  void closest_point(const int num_pts, const double *x, const double *y);
  std::vector<PolyApprox> get_polys();

private:
  int construct_tree(std::vector<KDCoord>::iterator pts_start, std::vector<KDCoord>::iterator pts_end, bool has_transformed, int level);
  double bb_sqr_dist(const int node_ind, const double x, const double y);
  double bb_sqr_dist(const MPIBB bb, const double x, const double y);
  bool check_if_locked_in(QueryPt &qp, const int num_ranks, const MPIBB *bb);
  void nearest_neighbour(double x, double y, int current_ind, std::vector<KDCoord>::iterator &closest_pt, double &closest_distance);

  std::set<int> cell_inds(std::vector<KDCoord> &points);
  void construct_polys(std::vector<KDCoord> &points, DGMesh *mesh, op_dat s);
  void update_poly_inds(std::vector<KDCoord> &points);

  std::vector<std::vector<KDCoord>::iterator> local_search(const int num_pts, const double *x, const double *y);
  std::vector<std::vector<KDCoord>::iterator> local_search(const int num_pts, const double *pts);

  void get_global_bounding_boxes(MPI_Comm *mpi_comm, MPIBB *mpi_bb);
  void round1_get_pts_to_send_to_ranks(const int num_pts, const double *x, const double *y,
                                       const MPIBB *mpi_bb, const int Reinit_comm_size,
                                       int *num_pts_to_send, std::map<int,std::vector<int>> &rank_to_pt_inds);
  void round1_prepare_send_recv(const int num_pts, const double *x, const double *y, 
                                const int Reinit_comm_size, std::map<int,std::vector<int>> &rank_to_pt_inds,
                                int *num_pts_to_recv, int *send_inds, int *recv_inds, 
                                double **pts_to_send, double **pts_to_recv, std::vector<int> &pt_send_rcv_map,
                                int &num_remote_pts);
  void round1_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm, 
                    int *num_pts_to_send, int *num_pts_to_recv, double *pts_to_send, double *pts_to_recv, 
                    int *send_inds, int *recv_inds, MPI_Request *requests);
  void round1_wait_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Request *requests, 
                         int *num_pts_to_send, int *num_pts_to_recv);
  void round1_send_results(const int Reinit_comm_rank, const int Reinit_comm_size, const int num_pts, const int num_remote_pts,
                           int *num_pts_to_send, int *num_pts_to_recv, int *send_inds, int *recv_inds,
                           std::vector<std::vector<KDCoord>::iterator> &remote_closest, MPIKDResponse **response,
                           MPI_Datatype *mpi_type, MPI_Comm *mpi_comm);
  std::vector<QueryPt> populate_query_pts(const int num_pts, const double *x, const double *y, const int Reinit_comm_rank, 
                                   const int Reinit_comm_size, MPIBB *mpi_bb, int *num_pts_to_send, MPIKDResponse *response, 
                                   std::vector<int> &pt_send_rcv_map, std::vector<std::vector<KDCoord>::iterator> &local_closest);
  void round2_pack_query_pts(const int Reinit_comm_size, int *num_pts_to_send, int *send_inds, 
                             std::map<int,std::vector<int>> &rank_to_qp, std::vector<QueryPt*> &nonLockedIn, 
                             double **round2_pts_to_send, std::vector<QueryPt*> &qp_ptrs);
  void round2_comms(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm, int *num_pts_to_send, int *num_pts_to_recv, 
                             int *send_inds, int *recv_inds, double *round2_pts_to_send, double **round2_pts_to_recv);
  void round2_results_comm(const int Reinit_comm_rank, const int Reinit_comm_size, MPI_Comm *mpi_comm, MPI_Datatype *mpi_type,
                           int *num_pts_to_send, int *num_pts_to_recv, int *send_inds, int *recv_inds, std::vector<std::vector<KDCoord>::iterator> &remote_closest,
                           MPIKDResponse **round2_send_response, MPIKDResponse **round2_recv_response);
  void round2_update_qp(const int Reinit_comm_size, int *num_pts_to_send, std::vector<QueryPt*> &qp_ptrs, MPIKDResponse *round2_recv_response);

  std::vector<KDNode> nodes;
  std::vector<KDCoord> points;
  int n;
  std::map<int,int> cell2polyMap;
  std::vector<PolyApprox> polys;
  static const int leaf_size = 10;
};

#endif
