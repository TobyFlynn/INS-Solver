#include "ls_utils/3d/kd_tree_mpi.h"

#include "timing.h"
#include "op2_utils.h"

#include "mpi.h"
#include "op_lib_mpi.h"

extern Timing *timer;

using namespace std;

KDTree3DMPI::KDTree3DMPI(DGMesh3D *m, const int alpha) : KDTree3D(m) {
  timer->startTimer("KDTree3DMPI - init");
  // Get local bounding box
  DG_FP min_x, max_x, min_y, max_y, min_z, max_z;
  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);
  const DG_FP *z_ptr = getOP2PtrHost(mesh->z, OP_READ);

  min_x = x_ptr[0];
  max_x = x_ptr[0];
  min_y = y_ptr[0];
  max_y = y_ptr[0];
  min_z = z_ptr[0];
  max_z = z_ptr[0];

  for(int i = 1; i < mesh->cells->size; i++) {
    if(x_ptr[i] < min_x) min_x = x_ptr[i];
    if(x_ptr[i] > min_x) max_x = x_ptr[i];
    if(y_ptr[i] < min_x) min_y = y_ptr[i];
    if(y_ptr[i] > min_x) max_y = y_ptr[i];
    if(z_ptr[i] < min_x) min_z = z_ptr[i];
    if(z_ptr[i] > min_x) max_z = z_ptr[i];
  }

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHost(mesh->z, OP_READ, z_ptr);

  // Get bounding boxes of all ranks
  int rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  DG_FP *bb_buf = (DG_FP *)calloc(comm_size * 6, sizeof(DG_FP));
  DG_FP snd_buf[] = {min_x, max_x, min_y, max_y, min_z, max_z};
  MPI_Allgather(snd_buf, 6, DG_MPI_FP, bb_buf, 6, DG_MPI_FP, MPI_COMM_WORLD);

  ranks.clear();
  const DG_FP min_l[] = {min_x, min_y, min_z};
  const DG_FP max_l[] = {max_x, max_y, max_z};
  for(int i = 0; i < comm_size; i++) {
    if(i == rank) continue;

    const DG_FP i_min_x = bb_buf[i * 6];
    const DG_FP i_max_x = bb_buf[i * 6 + 1];
    const DG_FP i_min_y = bb_buf[i * 6 + 2];
    const DG_FP i_max_y = bb_buf[i * 6 + 3];
    const DG_FP i_min_z = bb_buf[i * 6 + 4];
    const DG_FP i_max_z = bb_buf[i * 6 + 5];

    const DG_FP min_i[] = {bb_buf[i * 6], bb_buf[i * 6 + 2], bb_buf[i * 6 + 4]};
    const DG_FP max_i[] = {bb_buf[i * 6 + 1], bb_buf[i * 6 + 3], bb_buf[i * 6 + 5]};

    const DG_FP min_dist = min_dist_bb(min_l, max_l, min_i, max_i);

    if(alpha <= min_dist) {
      ranks.push_back(i);
    }
  }

  free(bb_buf);

  timer->endTimer("KDTree3DMPI - init");
}

DG_FP KDTree3DMPI::min_dist_bb(const DG_FP *min_0, const DG_FP *max_0,
                                    const DG_FP *min_1, const DG_FP *max_1) {
  DG_FP dists[] = {0.0, 0.0, 0.0};

  for(int d = 0; d < 3; d++) {
    if(!((min_1[d] >= min_0[d] && min_1[d] <= max_0[d]) || (max_1[d] >= min_0[d] && max_1[d] <= max_0[d]) || (min_1[d] <= min_0[d] && max_1[d] >= max_0[d]))) {
      if(min_1[d] > max_0[d]) {
        dists[d] = min_1[d] - max_0[d];
      } else if(max_1[d] < min_0[d]) {
        dists[d] = min_0[d] - max_1[d];
      } else {
        throw std::runtime_error("Logic error in KDTree min_dist_bb");
      }
    }
  }

  return sqrt(dists[0] * dists[0] + dists[1] * dists[1] + dists[2] * dists[2]);
}

void KDTree3DMPI::build_tree(const DG_FP *x, const DG_FP *y,
                                  const DG_FP *z, const int num, op_dat s) {
  pre_build_setup(x, y, z, num, s);

  timer->startTimer("KDTree3DMPI - comm");
  int rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

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
    MPI_Irecv(&recv_buffer[i * 2], 2, MPI_INT, ranks[i], 0, MPI_COMM_WORLD, &mpi_req_rcv[i]);
  }
  for(int i = 0; i < ranks.size(); i++) {
    MPI_Isend(num_pts_polys, 2, MPI_INT, ranks[i], 0, MPI_COMM_WORLD, &mpi_req_snd[i]);
  }
  MPI_Waitall(ranks.size(), mpi_req_rcv.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(ranks.size(), mpi_req_snd.data(), MPI_STATUS_IGNORE);

  int num_pts_g = 0;
  int num_polys_g = 0;
  for(int i = 0; i < ranks.size(); i++) {
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

  MPI_Waitall(ranks.size(), pts_mpi_req_rcv.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(ranks.size(), polys_mpi_req_rcv.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(ranks.size(), polys_coeff_mpi_req_rcv.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(ranks.size(), pts_mpi_req_snd.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(ranks.size(), polys_mpi_req_snd.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(ranks.size(), polys_coeff_mpi_req_snd.data(), MPI_STATUS_IGNORE);

  timer->endTimer("KDTree3DMPI - comm");

  timer->startTimer("KDTree3DMPI - unpack");
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
  timer->endTimer("KDTree3DMPI - unpack");

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
