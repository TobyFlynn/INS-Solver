#include "ls_utils/3d/kd_tree_mpi.h"

#include "timing.h"
#include "op2_utils.h"

#include "mpi.h"

extern Timing *timer;

using namespace std;

KDTree3DMPI::KDTree3DMPI(const DG_FP *x, const DG_FP *y, const DG_FP *z,
                         const int num, DGMesh3D *mesh, op_dat s) : KDTree3D(x, y, z, num, mesh, s) {

}

void KDTree3DMPI::build_tree() {
  int rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // Get number of points on each rank
  int num_pts_per_rank[comm_size];
  MPI_Allgather(&n, 1, MPI_INT, num_pts_per_rank, 1, MPI_INT, MPI_COMM_WORLD);
  int num_pts_g = 0;
  int num_pts_per_rank_disp[comm_size];
  for(int i = 0; i < comm_size; i++) {
    num_pts_per_rank_disp[i] = num_pts_g;
    num_pts_g += num_pts_per_rank[i];
  }

  // Get number of polynomial approximations per rank
  int num_polys_per_rank[comm_size];
  int num_polys = polys.size();
  MPI_Allgather(&num_polys, 1, MPI_INT, num_polys_per_rank, 1, MPI_INT, MPI_COMM_WORLD);
  int num_polys_g = 0;
  int num_polys_coeff_per_rank[comm_size];
  int num_polys_coeff_per_rank_disp[comm_size];
  for(int i = 0; i < comm_size; i++) {
    num_polys_coeff_per_rank_disp[i] = num_polys_g * PolyApprox3D::num_coeff();
    num_polys_g += num_polys_per_rank[i];
    num_polys_coeff_per_rank[i] = num_polys_per_rank[i] * PolyApprox3D::num_coeff();
  }

  // Prepare data to send
  vector<DG_FP> pts_x_snd(n);
  vector<DG_FP> pts_y_snd(n);
  vector<DG_FP> pts_z_snd(n);
  vector<int> pts_poly_snd(n);
  for(int i = 0; i < n; i++) {
    pts_x_snd[i] = points[i].x;
    pts_y_snd[i] = points[i].y;
    pts_z_snd[i] = points[i].z;
    pts_poly_snd[i] = points[i].poly;
  }

  vector<DG_FP> poly_coeff_snd(polys.size() * PolyApprox3D::num_coeff());
  for(int i = 0; i < polys.size(); i++) {
    int ind = i * PolyApprox3D::num_coeff();
    for(int j = 0; j < PolyApprox3D::num_coeff(); j++) {
      poly_coeff_snd[ind + j] = polys[i].get_coeff(j);
    }
  }

  // Send/Recieve
  vector<DG_FP> pts_x_rcv(num_pts_g);
  vector<DG_FP> pts_y_rcv(num_pts_g);
  vector<DG_FP> pts_z_rcv(num_pts_g);
  vector<int> pts_poly_rcv(num_pts_g);
  vector<DG_FP> poly_coeff_rcv(num_polys_g * PolyApprox3D::num_coeff());
  MPI_Allgatherv(pts_x_snd.data(), n, DG_MPI_FP, pts_x_rcv.data(), num_pts_per_rank, num_pts_per_rank_disp, DG_MPI_FP, MPI_COMM_WORLD);
  MPI_Allgatherv(pts_y_snd.data(), n, DG_MPI_FP, pts_y_rcv.data(), num_pts_per_rank, num_pts_per_rank_disp, DG_MPI_FP, MPI_COMM_WORLD);
  MPI_Allgatherv(pts_z_snd.data(), n, DG_MPI_FP, pts_z_rcv.data(), num_pts_per_rank, num_pts_per_rank_disp, DG_MPI_FP, MPI_COMM_WORLD);
  MPI_Allgatherv(pts_poly_snd.data(), n, MPI_INT, pts_poly_rcv.data(), num_pts_per_rank, num_pts_per_rank_disp, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgatherv(poly_coeff_snd.data(), num_polys * PolyApprox3D::num_coeff(), DG_MPI_FP, poly_coeff_rcv.data(), num_polys_coeff_per_rank, num_polys_coeff_per_rank_disp, DG_MPI_FP, MPI_COMM_WORLD);

  // Unpack data
  vector<KDCoord> newPts;
  int curr_rank = 0;
  for(int i = 0; i < num_pts_g; i++) {
    while(curr_rank < comm_size - 1 && i >= num_pts_per_rank_disp[curr_rank + 1])
      curr_rank++;
    KDCoord kc;
    kc.x = pts_x_rcv[i];
    kc.y = pts_y_rcv[i];
    kc.z = pts_z_rcv[i];
    kc.x_rot = pts_x_rcv[i];
    kc.y_rot = pts_y_rcv[i];
    kc.z_rot = pts_z_rcv[i];
    kc.poly = (num_polys_coeff_per_rank_disp[curr_rank] / PolyApprox3D::num_coeff()) + pts_poly_rcv[i];
    newPts.push_back(kc);
  }

  vector<PolyApprox3D> newPolys;
  for(int i = 0; i < num_polys_g; i++) {
    vector<DG_FP> c;
    int ind = i * PolyApprox3D::num_coeff();
    for(int j = 0; j < PolyApprox3D::num_coeff(); j++) {
      c.push_back(poly_coeff_rcv[ind + j]);
    }
    PolyApprox3D p(c, 0.0, 0.0, 0.0);
    newPolys.push_back(p);
  }

  points = newPts;
  polys = newPolys;

  if(points.size() == 0) {
    empty = true;
    return;
  }
  empty = false;

  timer->startTimer("K-D Tree - Construct Tree");
  construct_tree(points.begin(), points.end(), false, 0);
  timer->endTimer("K-D Tree - Construct Tree");
}