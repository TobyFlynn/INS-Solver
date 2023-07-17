#include "ls_utils/3d/kd_tree_mpi.h"

#include "timing.h"
#include "op2_utils.h"

#include "mpi.h"
#include "op_lib_mpi.h"

extern Timing *timer;

using namespace std;

KDTree3DMPI::KDTree3DMPI(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                         const int num, DGMesh3D *m, op_dat s) : KDTree3D(x, y, z, num, m, s) {

}

void KDTree3DMPI::build_tree() {
  int rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // Get ranks to communicate with
  halo_list exp_exec_list = OP_export_exec_list[mesh->cells->index];
  halo_list exp_nonexec_list = OP_export_nonexec_list[mesh->cells->index];
  halo_list imp_exec_list = OP_import_exec_list[mesh->cells->index];
  halo_list imp_nonexec_list = OP_import_nonexec_list[mesh->cells->index];
  set<int> ranks_set;
  for(int i = 0; i < exp_exec_list->ranks_size; i++) {
    ranks_set.insert(exp_exec_list->ranks[i]);
  }
  for(int i = 0; i < exp_nonexec_list->ranks_size; i++) {
    ranks_set.insert(exp_nonexec_list->ranks[i]);
  }
  for(int i = 0; i < imp_exec_list->ranks_size; i++) {
    ranks_set.insert(imp_exec_list->ranks[i]);
  }
  for(int i = 0; i < imp_nonexec_list->ranks_size; i++) {
    ranks_set.insert(imp_nonexec_list->ranks[i]);
  }
  vector<int> ranks;
  for(const int &r : ranks_set) {
    ranks.push_back(r);
  }

  // Get number of points and polynimials per rank
  vector<int> num_pts_per_rank(ranks.size());
  vector<int> num_polys_per_rank(ranks.size());
  vector<int> num_pts_per_rank_disp(ranks.size());
  vector<int> num_polys_per_rank_disp(ranks.size());
  vector<MPI_Request> mpi_req_snd(ranks.size());
  vector<MPI_Request> mpi_req_rcv(ranks.size());
  vector<int> recv_buffer(2 * ranks.size());
  int num_pts_polys[2] = {n, polys.size()};
  for(int i = 0; i < ranks.size(); i++) {
    MPI_Isend(num_pts_polys, 2, MPI_INT, ranks[i], 0, MPI_COMM_WORLD, &mpi_req_snd[i]);
  }
  for(int i = 0; i < ranks.size(); i++) {
    MPI_Irecv(&recv_buffer[i * 2], 2, MPI_INT, ranks[i], 0, MPI_COMM_WORLD, &mpi_req_rcv[i]);
  }
  int num_pts_g = 0;
  int num_polys_g = 0;
  for(int i = 0; i < ranks.size(); i++) {
    MPI_Wait(&mpi_req_rcv[i], MPI_STATUS_IGNORE);
    num_pts_per_rank[i] = recv_buffer[i * 2];
    num_polys_per_rank[i] = recv_buffer[i * 2 + 1];
    num_pts_per_rank_disp[i] = num_pts_g;
    num_polys_per_rank_disp[i] = num_polys_g;
    num_pts_g += num_pts_per_rank[i];
    num_polys_g += num_polys_per_rank[i];
  }

  // Prepare data to send
  vector<DG_FP> pts_snd(3 * n);
  vector<int> polys_snd(n);
  for(int i = 0; i < n; i++) {
    pts_snd[i * 3]     = points[i].x;
    pts_snd[i * 3 + 1] = points[i].y;
    pts_snd[i * 3 + 2] = points[i].z;
    polys_snd[i] = points[i].poly;
  }
  vector<DG_FP> poly_coeff_snd(polys.size() * PolyApprox3D::num_coeff());
  for(int i = 0; i < polys.size(); i++) {
    int ind = i * PolyApprox3D::num_coeff();
    for(int j = 0; j < PolyApprox3D::num_coeff(); j++) {
      poly_coeff_snd[ind + j] = polys[i].get_coeff(j);
    }
  }

  // Send/Recieve
  vector<DG_FP> pts_rcv(num_pts_g * 3);
  vector<int> polys_rcv(num_pts_g);
  vector<DG_FP> polys_coeff_rcv(num_polys_g * PolyApprox3D::num_coeff());
  vector<MPI_Request> pts_mpi_req_snd(ranks.size());
  vector<MPI_Request> pts_mpi_req_rcv(ranks.size());
  vector<MPI_Request> polys_mpi_req_snd(ranks.size());
  vector<MPI_Request> polys_mpi_req_rcv(ranks.size());
  vector<MPI_Request> polys_coeff_mpi_req_snd(ranks.size());
  vector<MPI_Request> polys_coeff_mpi_req_rcv(ranks.size());
  for(int i = 0; i < ranks.size(); i++) {
    MPI_Isend(pts_snd.data(), n * 3, DG_MPI_FP, ranks[i], 0, MPI_COMM_WORLD, &pts_mpi_req_snd[i]);
    MPI_Isend(polys_snd.data(), n, MPI_INT, ranks[i], 0, MPI_COMM_WORLD, &polys_mpi_req_snd[i]);
    MPI_Isend(poly_coeff_snd.data(), polys.size() * PolyApprox3D::num_coeff(), DG_MPI_FP, ranks[i], 0, MPI_COMM_WORLD, &polys_coeff_mpi_req_snd[i]);
  }
  for(int i = 0; i < ranks.size(); i++) {
    MPI_Irecv(&pts_rcv[num_pts_per_rank_disp[i] * 3], num_pts_per_rank[i] * 3, DG_MPI_FP, ranks[i], 0, MPI_COMM_WORLD, &pts_mpi_req_rcv[i]);
    MPI_Irecv(&polys_rcv[num_pts_per_rank_disp[i]], num_pts_per_rank[i], MPI_INT, ranks[i], 0, MPI_COMM_WORLD, &polys_mpi_req_rcv[i]);
    MPI_Irecv(&polys_coeff_rcv[num_polys_per_rank_disp[i] * PolyApprox3D::num_coeff()], num_polys_per_rank[i] * PolyApprox3D::num_coeff(), DG_MPI_FP, ranks[i], 0, MPI_COMM_WORLD, &polys_coeff_mpi_req_rcv[i]);
  }
  for(int i = 0; i < ranks.size(); i++) {
    MPI_Wait(&pts_mpi_req_rcv[i], MPI_STATUS_IGNORE);
    MPI_Wait(&polys_mpi_req_rcv[i], MPI_STATUS_IGNORE);
    MPI_Wait(&polys_coeff_mpi_req_rcv[i], MPI_STATUS_IGNORE);
  }

  // Unpack data
  int curr_rank = 0;
  int num_polys_l = polys.size();
  for(int i = 0; i < num_pts_g; i++) {
    while(curr_rank < ranks.size() - 1 && i >= num_pts_per_rank_disp[curr_rank + 1])
      curr_rank++;
    KDCoord kc;
    kc.x = pts_rcv[i * 3];
    kc.y = pts_rcv[i * 3 + 1];
    kc.z = pts_rcv[i * 3 + 2];
    kc.x_rot = pts_rcv[i * 3];
    kc.y_rot = pts_rcv[i * 3 + 1];
    kc.z_rot = pts_rcv[i * 3 + 2];
    kc.poly = num_polys_l + num_polys_per_rank_disp[curr_rank] + polys_rcv[i];
    points.push_back(kc);
  }
  for(int i = 0; i < num_polys_g; i++) {
    vector<DG_FP> c;
    int ind = i * PolyApprox3D::num_coeff();
    for(int j = 0; j < PolyApprox3D::num_coeff(); j++) {
      c.push_back(polys_coeff_rcv[ind + j]);
    }
    PolyApprox3D p(c, 0.0, 0.0, 0.0);
    polys.push_back(p);
  }

  if(points.size() == 0) {
    empty = true;
    // So saving timer info at end for MPI doesn't break
    timer->startTimer("K-D Tree - Construct Tree");
    timer->endTimer("K-D Tree - Construct Tree");
    return;
  }
  empty = false;

  timer->startTimer("K-D Tree - Construct Tree");
  construct_tree(points.begin(), points.end(), false, 0);
  timer->endTimer("K-D Tree - Construct Tree");
}
