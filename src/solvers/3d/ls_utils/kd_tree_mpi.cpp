#include "ls_utils/3d/kd_tree_mpi.h"

#include "timing.h"
#include "op2_utils.h"

#include "mpi.h"
#include "op_lib_mpi.h"

extern Timing *timer;

using namespace std;

struct kdtree_mpi_pack {
  DG_FP pos[3];
  int poly_ind;
};

KDTree3DMPI::KDTree3DMPI(DGMesh3D *m, const DG_FP alpha) : KDTree3D(m) {
  timer->startTimer("KDTree3DMPI - init");

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  const int block_lengths[] = {3,1};
  MPI_Aint displacements[2];
  MPI_Aint base_address;
  struct kdtree_mpi_pack dummy_packed;
  MPI_Get_address(&dummy_packed, &base_address);
  MPI_Get_address(&dummy_packed.pos[0], &displacements[0]);
  MPI_Get_address(&dummy_packed.poly_ind, &displacements[1]);
  displacements[0] = MPI_Aint_diff(displacements[0], base_address);
  displacements[1] = MPI_Aint_diff(displacements[1], base_address);
  MPI_Datatype mpi_datatypes[] = {DG_MPI_FP, MPI_INT};
  MPI_Type_create_struct(2, block_lengths, displacements, mpi_datatypes, &packed_type);

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

  for(int i = 1; i < mesh->cells->size * DG_NP; i++) {
    if(x_ptr[i] < min_x) min_x = x_ptr[i];
    if(x_ptr[i] > max_x) max_x = x_ptr[i];
    if(y_ptr[i] < min_x) min_y = y_ptr[i];
    if(y_ptr[i] > max_y) max_y = y_ptr[i];
    if(z_ptr[i] < min_x) min_z = z_ptr[i];
    if(z_ptr[i] > max_z) max_z = z_ptr[i];
  }

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHost(mesh->z, OP_READ, z_ptr);

  // Get bounding boxes of all ranks
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

    const DG_FP min_i[] = {i_min_x, i_min_y, i_min_z};
    const DG_FP max_i[] = {i_max_x, i_max_y, i_max_z};

    const DG_FP min_dist = min_dist_bb(min_l, max_l, min_i, max_i);

    if(min_dist <= alpha) {
      ranks.push_back(i);
    }
  }

  free(bb_buf);

  timer->endTimer("KDTree3DMPI - init");
}

DG_FP KDTree3DMPI::min_dist_bb(const DG_FP *min_0, const DG_FP *max_0,
                               const DG_FP *min_1, const DG_FP *max_1) {
/*
  DG_FP a[] = {
                fmax(0.0, min_0[0] - max_1[0]),
                fmax(0.0, min_0[1] - max_1[1]),
                fmax(0.0, min_0[2] - max_1[2])
              };
  DG_FP b[] = {
                fmax(0.0, min_1[0] - max_0[0]),
                fmax(0.0, min_1[1] - max_0[1]),
                fmax(0.0, min_1[2] - max_0[2])
              };

  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
*/

  DG_FP min_i[] = {
    fmax(min_0[0], min_1[0]),
    fmax(min_0[1], min_1[1]),
    fmax(min_0[2], min_1[2])
  };

  DG_FP max_i[] = {
    fmin(max_0[0], max_1[0]),
    fmin(max_0[1], max_1[1]),
    fmin(max_0[2], max_1[2])
  };

  DG_FP dist = 0.0;
  for(int i = 0; i < 3; i++) {
    if(min_i[i] > max_i[i])
      dist += (min_i[i] - max_i[i]) * (min_i[i] - max_i[i]);
  }
  return sqrt(dist);
}

void KDTree3DMPI::build_tree(const DG_FP *x, const DG_FP *y,
                                  const DG_FP *z, const int num, op_dat s) {
  timer->startTimer("KDTree3DMPI - pre_build_setup");
  pre_build_setup(x, y, z, num, s);
  timer->endTimer("KDTree3DMPI - pre_build_setup");

  timer->startTimer("KDTree3DMPI - comm - 0");
  // Get number of points and polynimials per rank
  vector<int> num_pts_per_rank(ranks.size());
  vector<int> num_polys_per_rank(ranks.size());
  vector<int> num_pts_per_rank_disp(ranks.size());
  vector<int> num_polys_per_rank_disp(ranks.size());
  vector<MPI_Request> mpi_req_snd(ranks.size());
  vector<MPI_Request> mpi_req_rcv(ranks.size());
  vector<int> recv_buffer(2 * ranks.size());
  int num_pts_polys[2] = {n, static_cast<int>(polys.size())};
  for(int i = 0; i < ranks.size(); i++) {
    MPI_Irecv(&recv_buffer[i * 2], 2, MPI_INT, ranks[i], 0, MPI_COMM_WORLD, &mpi_req_rcv[i]);
  }
  for(int i = 0; i < ranks.size(); i++) {
    MPI_Isend(num_pts_polys, 2, MPI_INT, ranks[i], 0, MPI_COMM_WORLD, &mpi_req_snd[i]);
  }
  MPI_Waitall(ranks.size(), mpi_req_rcv.data(), MPI_STATUS_IGNORE);

  timer->endTimer("KDTree3DMPI - comm - 0");
  timer->startTimer("KDTree3DMPI - packing");

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
  vector<kdtree_mpi_pack> packed_data(n);
  for(int i = 0; i < n; i++) {
    packed_data[i].pos[0] = points[i].x;
    packed_data[i].pos[1] = points[i].y;
    packed_data[i].pos[2] = points[i].z;
    packed_data[i].poly_ind = points[i].poly;
  }
  vector<DG_FP> poly_coeff_snd(polys.size() * PolyApprox3D::num_coeff());
  for(int i = 0; i < polys.size(); i++) {
    int ind = i * PolyApprox3D::num_coeff();
    for(int j = 0; j < PolyApprox3D::num_coeff(); j++) {
      poly_coeff_snd[ind + j] = polys[i].get_coeff(j);
    }
  }
  timer->endTimer("KDTree3DMPI - packing");

  timer->startTimer("KDTree3DMPI - comm - 1");
  // Send/Recieve
  vector<kdtree_mpi_pack> pts_rcv(num_pts_g);
  vector<DG_FP> polys_coeff_rcv(num_polys_g * PolyApprox3D::num_coeff());
  vector<MPI_Request> pts_mpi_req_snd(ranks.size());
  vector<MPI_Request> pts_mpi_req_rcv(ranks.size());
  vector<MPI_Request> polys_coeff_mpi_req_snd(ranks.size());
  vector<MPI_Request> polys_coeff_mpi_req_rcv(ranks.size());

  for(int i = 0; i < ranks.size(); i++) {
    MPI_Irecv(&pts_rcv[num_pts_per_rank_disp[i]], num_pts_per_rank[i], packed_type, ranks[i], 0, MPI_COMM_WORLD, &pts_mpi_req_rcv[i]);
    MPI_Irecv(&polys_coeff_rcv[num_polys_per_rank_disp[i] * PolyApprox3D::num_coeff()], num_polys_per_rank[i] * PolyApprox3D::num_coeff(), DG_MPI_FP, ranks[i], 0, MPI_COMM_WORLD, &polys_coeff_mpi_req_rcv[i]);
  }

  for(int i = 0; i < ranks.size(); i++) {
    MPI_Isend(packed_data.data(), n, packed_type, ranks[i], 0, MPI_COMM_WORLD, &pts_mpi_req_snd[i]);
    MPI_Isend(poly_coeff_snd.data(), polys.size() * PolyApprox3D::num_coeff(), DG_MPI_FP, ranks[i], 0, MPI_COMM_WORLD, &polys_coeff_mpi_req_snd[i]);
  }

  MPI_Waitall(ranks.size(), pts_mpi_req_rcv.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(ranks.size(), polys_coeff_mpi_req_rcv.data(), MPI_STATUS_IGNORE);

  timer->endTimer("KDTree3DMPI - comm - 1");

  timer->startTimer("KDTree3DMPI - unpack");
  // Unpack data
  int curr_rank = 0;
  int num_polys_l = polys.size();
  for(int i = 0; i < num_pts_g; i++) {
    while(curr_rank < ranks.size() - 1 && i >= num_pts_per_rank_disp[curr_rank + 1])
      curr_rank++;
    KDCoord kc;
    kc.x = pts_rcv[i].pos[0];
    kc.y = pts_rcv[i].pos[1];
    kc.z = pts_rcv[i].pos[2];
    kc.x_rot = pts_rcv[i].pos[0];
    kc.y_rot = pts_rcv[i].pos[1];
    kc.z_rot = pts_rcv[i].pos[2];
    kc.poly = num_polys_l + num_polys_per_rank_disp[curr_rank] + pts_rcv[i].poly_ind;
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
  } else {
    empty = false;

    timer->startTimer("K-D Tree - Construct Tree");
    construct_tree(points.begin(), points.end(), false, 0);
    timer->endTimer("K-D Tree - Construct Tree");
  }

  timer->startTimer("KDTree3DMPI - comm - wait snd");
  MPI_Waitall(ranks.size(), mpi_req_snd.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(ranks.size(), pts_mpi_req_snd.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(ranks.size(), polys_coeff_mpi_req_snd.data(), MPI_STATUS_IGNORE);
  timer->endTimer("KDTree3DMPI - comm - wait snd");

  int min_n, max_n, min_ranks;
  int min_pol, max_pol, max_ranks;
  int num_ranks = ranks.size();
  MPI_Reduce(&num_pts_polys[0], &min_n, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&num_pts_polys[0], &max_n, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&num_pts_polys[1], &min_pol, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&num_pts_polys[1], &max_pol, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&num_ranks, &min_ranks, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&num_ranks, &max_ranks, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  // op_printf("min n: %d max n: %d\n", min_n, max_n);
  // op_printf("min pol: %d max pol: %d\n", min_pol, max_pol);
  // op_printf("min ranks: %d max ranks: %d\n", min_ranks, max_ranks);
}
