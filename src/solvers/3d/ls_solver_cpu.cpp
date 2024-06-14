#include "solvers/3d/ls_solver.h"

#include <set>
#include <random>

#include "dg_op2_blas.h"
#include "dg_utils.h"
#include "op2_utils.h"
#include "ls_utils/3d/poly_approx.h"
#include "dg_dat_pool.h"
#include "dg_constants/dg_constants_3d.h"
#include "dg_global_constants/dg_global_constants_3d.h"

extern DGDatPool *dg_dat_pool;
extern DGConstants *constants;

void rst2xyz(DG_FP &sampleX, DG_FP &sampleY, DG_FP &sampleZ,
             const DG_FP *nodeX, const DG_FP *nodeY,
             const DG_FP *nodeZ) {
  DG_FP r_ = sampleX;
  DG_FP s_ = sampleY;
  DG_FP t_ = sampleZ;

  sampleX = 0.5 * (-(1.0 + r_ + s_ + t_) * nodeX[0] + (1.0 + r_) * nodeX[1] + (1.0 + s_) * nodeX[2] + (1.0 + t_) * nodeX[3]);
  sampleY = 0.5 * (-(1.0 + r_ + s_ + t_) * nodeY[0] + (1.0 + r_) * nodeY[1] + (1.0 + s_) * nodeY[2] + (1.0 + t_) * nodeY[3]);
  sampleZ = 0.5 * (-(1.0 + r_ + s_ + t_) * nodeZ[0] + (1.0 + r_) * nodeZ[1] + (1.0 + s_) * nodeZ[2] + (1.0 + t_) * nodeZ[3]);
}

bool pt_in_tetra(const DG_FP ptX, const DG_FP ptY, const DG_FP ptZ,
                 const DG_FP *nodeX, const DG_FP *nodeY,
                 const DG_FP *nodeZ) {
  bool sameSide0, sameSide1, sameSide2, sameSide3;
  DG_FP normal[3];
  // (v1 - v0) x (v2 - v0)
  normal[0] = (nodeY[1] - nodeY[0]) * (nodeZ[2] - nodeZ[0]) - (nodeZ[1] - nodeZ[0]) * (nodeY[2] - nodeY[0]);
  normal[1] = (nodeZ[1] - nodeZ[0]) * (nodeX[2] - nodeX[0]) - (nodeX[1] - nodeX[0]) * (nodeZ[2] - nodeZ[0]);
  normal[2] = (nodeX[1] - nodeX[0]) * (nodeY[2] - nodeY[0]) - (nodeY[1] - nodeY[0]) * (nodeX[2] - nodeX[0]);
  // normal . (v3 - v0)
  DG_FP dotV = normal[0] * (nodeX[3] - nodeX[0]) + normal[1] * (nodeY[3] - nodeY[0]) + normal[2] * (nodeZ[3] - nodeZ[0]);
  // normal . (p - v0)
  DG_FP dotP = normal[0] * (ptX - nodeX[0]) + normal[1] * (ptY - nodeY[0]) + normal[2] * (ptZ - nodeZ[0]);
  sameSide0 = (dotV > 0.0) == (dotP > 0.0);

  // (v2 - v1) x (v3 - v1)
  normal[0] = (nodeY[2] - nodeY[1]) * (nodeZ[3] - nodeZ[1]) - (nodeZ[2] - nodeZ[1]) * (nodeY[3] - nodeY[1]);
  normal[1] = (nodeZ[2] - nodeZ[1]) * (nodeX[3] - nodeX[1]) - (nodeX[2] - nodeX[1]) * (nodeZ[3] - nodeZ[1]);
  normal[2] = (nodeX[2] - nodeX[1]) * (nodeY[3] - nodeY[1]) - (nodeY[2] - nodeY[1]) * (nodeX[3] - nodeX[1]);
  // normal . (v0 - v1)
  dotV = normal[0] * (nodeX[0] - nodeX[1]) + normal[1] * (nodeY[0] - nodeY[1]) + normal[2] * (nodeZ[0] - nodeZ[1]);
  // normal . (p - v1)
  dotP = normal[0] * (ptX - nodeX[1]) + normal[1] * (ptY - nodeY[1]) + normal[2] * (ptZ - nodeZ[1]);
  sameSide1 = (dotV > 0.0) == (dotP > 0.0);

  // (v3 - v2) x (v0 - v2)
  normal[0] = (nodeY[3] - nodeY[2]) * (nodeZ[0] - nodeZ[2]) - (nodeZ[3] - nodeZ[2]) * (nodeY[0] - nodeY[2]);
  normal[1] = (nodeZ[3] - nodeZ[2]) * (nodeX[0] - nodeX[2]) - (nodeX[3] - nodeX[2]) * (nodeZ[0] - nodeZ[2]);
  normal[2] = (nodeX[3] - nodeX[2]) * (nodeY[0] - nodeY[2]) - (nodeY[3] - nodeY[2]) * (nodeX[0] - nodeX[2]);
  // normal . (v1 - v2)
  dotV = normal[0] * (nodeX[1] - nodeX[2]) + normal[1] * (nodeY[1] - nodeY[2]) + normal[2] * (nodeZ[1] - nodeZ[2]);
  // normal . (p - v2)
  dotP = normal[0] * (ptX - nodeX[2]) + normal[1] * (ptY - nodeY[2]) + normal[2] * (ptZ - nodeZ[2]);
  sameSide2 = (dotV > 0.0) == (dotP > 0.0);

  // (v0 - v3) x (v1 - v3)
  normal[0] = (nodeY[0] - nodeY[3]) * (nodeZ[1] - nodeZ[3]) - (nodeZ[0] - nodeZ[3]) * (nodeY[1] - nodeY[3]);
  normal[1] = (nodeZ[0] - nodeZ[3]) * (nodeX[1] - nodeX[3]) - (nodeX[0] - nodeX[3]) * (nodeZ[1] - nodeZ[3]);
  normal[2] = (nodeX[0] - nodeX[3]) * (nodeY[1] - nodeY[3]) - (nodeY[0] - nodeY[3]) * (nodeX[1] - nodeX[3]);
  // normal . (v2 - v3)
  dotV = normal[0] * (nodeX[2] - nodeX[3]) + normal[1] * (nodeY[2] - nodeY[3]) + normal[2] * (nodeZ[2] - nodeZ[3]);
  // normal . (p - v3)
  dotP = normal[0] * (ptX - nodeX[3]) + normal[1] * (ptY - nodeY[3]) + normal[2] * (ptZ - nodeZ[3]);
  sameSide3 = (dotV > 0.0) == (dotP > 0.0);

  return sameSide0 && sameSide1 && sameSide2 && sameSide3;
}

template<int NUM_PER_LINE>
void set_sample_start_coords(DG_FP *r, DG_FP *s, DG_FP *t) {
  DG_FP node0[] = {-1, -1, -1};
  DG_FP node1[] = {1, -1, -1};
  DG_FP node2[] = {-1, 1, -1};
  DG_FP node3[] = {-1, -1, 1};

  DG_FP r_c = (node0[0] + node1[0] + node2[0] + node3[0]) / 4.0;
  DG_FP s_c = (node0[1] + node1[1] + node2[1] + node3[1]) / 4.0;
  DG_FP t_c = (node0[2] + node1[2] + node2[2] + node3[2]) / 4.0;

  r[0] = r_c;
  s[0] = s_c;
  t[0] = t_c;

  // Line from centre to node0
  for(int i = 1; i < NUM_PER_LINE; i++) {
    const DG_FP frac = (DG_FP)i / (DG_FP)NUM_PER_LINE;
    r[i] = frac * (node0[0] - r_c) + r_c;
    s[i] = frac * (node0[1] - s_c) + s_c;
    t[i] = frac * (node0[2] - t_c) + t_c;
  }
  r[NUM_PER_LINE] = node0[0];
  s[NUM_PER_LINE] = node0[1];
  t[NUM_PER_LINE] = node0[2];

  // Line from centre to node1
  for(int i = 1; i < NUM_PER_LINE; i++) {
    const DG_FP frac = (DG_FP)i / (DG_FP)NUM_PER_LINE;
    r[NUM_PER_LINE + i] = frac * (node1[0] - r_c) + r_c;
    s[NUM_PER_LINE + i] = frac * (node1[1] - s_c) + s_c;
    t[NUM_PER_LINE + i] = frac * (node1[2] - t_c) + t_c;
  }
  r[2 * NUM_PER_LINE] = node1[0];
  s[2 * NUM_PER_LINE] = node1[1];
  t[2 * NUM_PER_LINE] = node1[2];

  // Line from centre to node2
  for(int i = 1; i < NUM_PER_LINE; i++) {
    const DG_FP frac = (DG_FP)i / (DG_FP)NUM_PER_LINE;
    r[2 * NUM_PER_LINE + i] = frac * (node2[0] - r_c) + r_c;
    s[2 * NUM_PER_LINE + i] = frac * (node2[1] - s_c) + s_c;
    t[2 * NUM_PER_LINE + i] = frac * (node2[2] - t_c) + t_c;
  }
  r[3 * NUM_PER_LINE] = node2[0];
  s[3 * NUM_PER_LINE] = node2[1];
  t[3 * NUM_PER_LINE] = node2[2];

  // Line from centre to node3
  for(int i = 1; i < NUM_PER_LINE; i++) {
    const DG_FP frac = (DG_FP)i / (DG_FP)NUM_PER_LINE;
    r[3 * NUM_PER_LINE + i] = frac * (node3[0] - r_c) + r_c;
    s[3 * NUM_PER_LINE + i] = frac * (node3[1] - s_c) + s_c;
    t[3 * NUM_PER_LINE + i] = frac * (node3[2] - t_c) + t_c;
  }
  r[4 * NUM_PER_LINE] = node3[0];
  s[4 * NUM_PER_LINE] = node3[1];
  t[4 * NUM_PER_LINE] = node3[2];
}

void intersecting_pts(const DG_FP *s, const DG_FP *nodeX, const DG_FP *nodeY,
                      const DG_FP *nodeZ, std::vector<DGUtils::Vec<3>> &pts) {
  // Check each tet edge
  // TODO generalise for all orders
  const int fmask_node_ind_0 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int fmask_node_ind_1 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + DG_NPF - 1];
  const int fmask_node_ind_2 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + 2 * DG_NPF - 1];
  const int fmask_node_ind_3 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + 3 * DG_NPF - 1];

  const DG_FP node0_s = s[0];
  const DG_FP node1_s = s[3];
  const DG_FP node2_s = s[9];
  const DG_FP node3_s = s[19];
  // Node0 -> Node1
  if(node0_s > 0.0 != node1_s > 0.0) {
    const DG_FP s_factor = node0_s / (node0_s - node1_s);
    DGUtils::Vec<3> coord;
    coord[0] = nodeX[0] - s_factor * (nodeX[0] - nodeX[1]);
    coord[1] = nodeY[0] - s_factor * (nodeY[0] - nodeY[1]);
    coord[2] = nodeZ[0] - s_factor * (nodeZ[0] - nodeZ[1]);
    pts.push_back(coord);
  }
  // Node0 -> Node2
  if(node0_s > 0.0 != node2_s > 0.0) {
    const DG_FP s_factor = node0_s / (node0_s - node2_s);
    DGUtils::Vec<3> coord;
    coord[0] = nodeX[0] - s_factor * (nodeX[0] - nodeX[2]);
    coord[1] = nodeY[0] - s_factor * (nodeY[0] - nodeY[2]);
    coord[2] = nodeZ[0] - s_factor * (nodeZ[0] - nodeZ[2]);
    pts.push_back(coord);
  }
  // Node0 -> Node3
  if(node0_s > 0.0 != node3_s > 0.0) {
    const DG_FP s_factor = node0_s / (node0_s - node3_s);
    DGUtils::Vec<3> coord;
    coord[0] = nodeX[0] - s_factor * (nodeX[0] - nodeX[3]);
    coord[1] = nodeY[0] - s_factor * (nodeY[0] - nodeY[3]);
    coord[2] = nodeZ[0] - s_factor * (nodeZ[0] - nodeZ[3]);
    pts.push_back(coord);
  }
  // Node1 -> Node2
  if(node1_s > 0.0 != node2_s > 0.0) {
    const DG_FP s_factor = node1_s / (node1_s - node2_s);
    DGUtils::Vec<3> coord;
    coord[0] = nodeX[1] - s_factor * (nodeX[1] - nodeX[2]);
    coord[1] = nodeY[1] - s_factor * (nodeY[1] - nodeY[2]);
    coord[2] = nodeZ[1] - s_factor * (nodeZ[1] - nodeZ[2]);
    pts.push_back(coord);
  }
  // Node1 -> Node3
  if(node1_s > 0.0 != node3_s > 0.0) {
    const DG_FP s_factor = node1_s / (node1_s - node3_s);
    DGUtils::Vec<3> coord;
    coord[0] = nodeX[1] - s_factor * (nodeX[1] - nodeX[3]);
    coord[1] = nodeY[1] - s_factor * (nodeY[1] - nodeY[3]);
    coord[2] = nodeZ[1] - s_factor * (nodeZ[1] - nodeZ[3]);
    pts.push_back(coord);
  }
  // Node2 -> Node3
  if(node2_s > 0.0 != node3_s > 0.0) {
    const DG_FP s_factor = node2_s / (node2_s - node3_s);
    DGUtils::Vec<3> coord;
    coord[0] = nodeX[2] - s_factor * (nodeX[2] - nodeX[3]);
    coord[1] = nodeY[2] - s_factor * (nodeY[2] - nodeY[3]);
    coord[2] = nodeZ[2] - s_factor * (nodeZ[2] - nodeZ[3]);
    pts.push_back(coord);
  }
}

void intersect_3pts(const std::vector<DGUtils::Vec<3>> &intersect_pts, DG_FP *sampleX,
                    DG_FP *sampleY, DG_FP *sampleZ) {
  // Get intersecting plane eqn
  DG_FP plane_coeff[4];
  DGUtils::Vec<3> vec0 = intersect_pts[1] - intersect_pts[0];
  DGUtils::Vec<3> vec1 = intersect_pts[2] - intersect_pts[0];

  plane_coeff[0] = vec0[1] * vec1[2] - vec0[2] * vec1[1];
  plane_coeff[1] = vec0[2] * vec1[0] - vec0[0] * vec1[2];
  plane_coeff[2] = vec0[0] * vec1[1] - vec0[1] * vec1[0];
  plane_coeff[3] = -(plane_coeff[0] * intersect_pts[0][0] + plane_coeff[1] * intersect_pts[0][1] + plane_coeff[2] * intersect_pts[0][2]);

  // Calc sampling points
  DG_FP r_[] = {-1, -0.333333333333333, 0.333333333333333, 1, -1, -0.333333333333333, 0.333333333333333, -1, -0.333333333333333, -1};
  DG_FP s_[] = {-1, -1, -1, -1, -0.333333333333333, -0.333333333333333, -0.333333333333333, 0.333333333333333, 0.333333333333333, 1};

  for(int i = 0; i < 10; i++) {
    r_[i] *= 0.75;
    s_[i] *= 0.75;
  }

  DG_FP x_2D[10], y_2D[10], z_2D[10];
  for(int i = 0; i < 10; i++) {
    x_2D[i] = 0.5*(-(r_[i]+s_[i])*intersect_pts[0][0]+(1+r_[i])*intersect_pts[1][0]+(1+s_[i])*intersect_pts[2][0]);
    y_2D[i] = 0.5*(-(r_[i]+s_[i])*intersect_pts[0][1]+(1+r_[i])*intersect_pts[1][1]+(1+s_[i])*intersect_pts[2][1]);
    z_2D[i] = 0.5*(-(r_[i]+s_[i])*intersect_pts[0][2]+(1+r_[i])*intersect_pts[1][2]+(1+s_[i])*intersect_pts[2][2]);
  }

  // Now project onto intersecting 3D plane
  for(int i = 0; i < 10; i++) {
    if(fabs(plane_coeff[2]) > 1e-8) {
      sampleX[i] = x_2D[i];
      sampleY[i] = y_2D[i];
      sampleZ[i] = (-plane_coeff[0] * x_2D[i] - plane_coeff[1] * y_2D[i] - plane_coeff[3]) / plane_coeff[2];
    } else if(fabs(plane_coeff[1]) > 1e-8) {
      sampleX[i] = x_2D[i];
      sampleY[i] = (-plane_coeff[0] * x_2D[i] - plane_coeff[2] * z_2D[i] - plane_coeff[3]) / plane_coeff[1];
      sampleZ[i] = z_2D[i];
    } else if(fabs(plane_coeff[0]) > 1e-8) {
      sampleX[i] = (-plane_coeff[1] * y_2D[i] - plane_coeff[2] * z_2D[i] - plane_coeff[3]) / plane_coeff[0];
      sampleY[i] = y_2D[i];
      sampleZ[i] = z_2D[i];
    } else {
      sampleX[i] = NAN;
      sampleY[i] = NAN;
      sampleZ[i] = NAN;
    }
  }
}

void intersect_4pts(const std::vector<DGUtils::Vec<3>> &intersect_pts, DG_FP *sampleX,
                    DG_FP *sampleY, DG_FP *sampleZ) {
  // 1st Triangle is first 3 points
  intersect_3pts(intersect_pts, sampleX, sampleY, sampleZ);
  // Find which points make up second triangle
  DGUtils::Vec<3> mid_0_1 = 0.5 * (intersect_pts[0] + intersect_pts[1]);
  DGUtils::Vec<3> mid_1_2 = 0.5 * (intersect_pts[1] + intersect_pts[2]);
  DGUtils::Vec<3> mid_2_0 = 0.5 * (intersect_pts[2] + intersect_pts[0]);

  DG_FP dist_0_1 = (mid_0_1 - intersect_pts[3]).sqr_magnitude();
  DG_FP dist_1_2 = (mid_1_2 - intersect_pts[3]).sqr_magnitude();
  DG_FP dist_2_0 = (mid_2_0 - intersect_pts[3]).sqr_magnitude();

  // Call 2nd Triangle
  if(dist_0_1 <= dist_1_2 && dist_0_1 <= dist_2_0) {
    std::vector<DGUtils::Vec<3>> tmp_intersect_pts;
    tmp_intersect_pts.push_back(intersect_pts[0]);
    tmp_intersect_pts.push_back(intersect_pts[1]);
    tmp_intersect_pts.push_back(intersect_pts[3]);
    intersect_3pts(tmp_intersect_pts, sampleX + 10, sampleY + 10, sampleZ + 10);
  } else if(dist_1_2 <= dist_0_1 && dist_1_2 <= dist_2_0) {
    std::vector<DGUtils::Vec<3>> tmp_intersect_pts;
    tmp_intersect_pts.push_back(intersect_pts[1]);
    tmp_intersect_pts.push_back(intersect_pts[2]);
    tmp_intersect_pts.push_back(intersect_pts[3]);
    intersect_3pts(tmp_intersect_pts, sampleX + 10, sampleY + 10, sampleZ + 10);
  } else {
    std::vector<DGUtils::Vec<3>> tmp_intersect_pts;
    tmp_intersect_pts.push_back(intersect_pts[2]);
    tmp_intersect_pts.push_back(intersect_pts[0]);
    tmp_intersect_pts.push_back(intersect_pts[3]);
    intersect_3pts(tmp_intersect_pts, sampleX + 10, sampleY + 10, sampleZ + 10);
  }
}

bool simplified_newton(DG_FP &pt_x, DG_FP &pt_y, DG_FP &pt_z, PolyApprox3D &pol, const DG_FP tol) {
  bool converged = false;
  for(int step = 0; step < 100; step++) {
    DG_FP surf = pol.val_at(pt_x, pt_y, pt_z);
    DG_FP dsdx, dsdy, dsdz;
    pol.grad_at(pt_x, pt_y, pt_z, dsdx, dsdy, dsdz);

    DG_FP sqrnorm = dsdx * dsdx + dsdy * dsdy + dsdz * dsdz;
    if(sqrnorm > 1e-14) {
      dsdx *= surf / sqrnorm;
      dsdy *= surf / sqrnorm;
      dsdz *= surf / sqrnorm;
    }

    pt_x -= dsdx;
    pt_y -= dsdy;
    pt_z -= dsdz;

    // Check convergence
    if(dsdx * dsdx + dsdy * dsdy + dsdz * dsdz < tol) {
      converged = true;
      break;
    }
  }
  return converged;
}

void LevelSetSolver3D::sampleInterface(op_dat sampleX, op_dat sampleY, op_dat sampleZ,
              std::vector<PolyApprox3D> &polys, std::map<int,int> &cell2polyMap,
              std::set<int> &cellInds) {
  // Setup random number generator for later
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1.0, 1.0);

  DG_FP ref_r_ptr[LS_SAMPLE_NP], ref_s_ptr[LS_SAMPLE_NP], ref_t_ptr[LS_SAMPLE_NP];
  // set_sample_start_coords<10>(ref_r_ptr, ref_s_ptr, ref_t_ptr);
  const DG_FP *tmp_r_ptr = constants->get_mat_ptr(DGConstants::R) + (DG_ORDER - 1) * DG_NP;
  const DG_FP *tmp_s_ptr = constants->get_mat_ptr(DGConstants::S) + (DG_ORDER - 1) * DG_NP;
  const DG_FP *tmp_t_ptr = constants->get_mat_ptr(DGConstants::T) + (DG_ORDER - 1) * DG_NP;

  for(int i = 0; i < DG_NP; i++) {
    ref_r_ptr[i] = tmp_r_ptr[i] * 0.75;
    ref_s_ptr[i] = tmp_s_ptr[i] * 0.75;
    ref_t_ptr[i] = tmp_t_ptr[i] * 0.75;
  }

  // const DG_FP *ref_r_ptr = constants->get_mat_ptr(DGConstants::R) + (DG_ORDER - 1) * DG_NP;
  // const DG_FP *ref_s_ptr = constants->get_mat_ptr(DGConstants::S) + (DG_ORDER - 1) * DG_NP;
  // const DG_FP *ref_t_ptr = constants->get_mat_ptr(DGConstants::T) + (DG_ORDER - 1) * DG_NP;

  const DG_FP *s_ptr = getOP2PtrHost(s, OP_READ);
  const DG_FP *nodeX_ptr = getOP2PtrHost(mesh->nodeX, OP_READ);
  const DG_FP *nodeY_ptr = getOP2PtrHost(mesh->nodeY, OP_READ);
  const DG_FP *nodeZ_ptr = getOP2PtrHost(mesh->nodeZ, OP_READ);
  const DG_FP *x_ptr = getOP2PtrHost(mesh->x, OP_READ);
  const DG_FP *y_ptr = getOP2PtrHost(mesh->y, OP_READ);
  const DG_FP *z_ptr = getOP2PtrHost(mesh->z, OP_READ);
  DG_FP *sampleX_ptr = getOP2PtrHost(sampleX, OP_WRITE);
  DG_FP *sampleY_ptr = getOP2PtrHost(sampleY, OP_WRITE);
  DG_FP *sampleZ_ptr = getOP2PtrHost(sampleZ, OP_WRITE);

  #pragma omp parallel for
  for(int cell = 0; cell < mesh->cells->size; cell++) {
    DG_FP *sampleX_c = sampleX_ptr + cell * LS_SAMPLE_NP;
    DG_FP *sampleY_c = sampleY_ptr + cell * LS_SAMPLE_NP;
    DG_FP *sampleZ_c = sampleZ_ptr + cell * LS_SAMPLE_NP;
    for(int i = 0; i < LS_SAMPLE_NP; i++) {
      sampleX_c[i] = NAN;
      sampleY_c[i] = NAN;
      sampleZ_c[i] = NAN;
    }
  }

  std::vector<int> cell_inds_vec(cellInds.begin(), cellInds.end());
  const DG_FP tol = fmax(1e-18, 1e-4 * h * h);
  // const DG_FP tol = 1e-16;
  // printf("TOL: %g\n", tol);

  #pragma omp parallel for
  for(int i = 0; i < cell_inds_vec.size(); i++) {
    const int cell = cell_inds_vec[i];
    const DG_FP *s_c = s_ptr + cell * DG_NP;
    const DG_FP *nodeX_c = nodeX_ptr + cell * 4;
    const DG_FP *nodeY_c = nodeY_ptr + cell * 4;
    const DG_FP *nodeZ_c = nodeZ_ptr + cell * 4;
    const DG_FP *x_c = x_ptr + cell * DG_NP;
    const DG_FP *y_c = y_ptr + cell * DG_NP;
    const DG_FP *z_c = z_ptr + cell * DG_NP;
    DG_FP *sampleX_c = sampleX_ptr + cell * LS_SAMPLE_NP;
    DG_FP *sampleY_c = sampleY_ptr + cell * LS_SAMPLE_NP;
    DG_FP *sampleZ_c = sampleZ_ptr + cell * LS_SAMPLE_NP;

    // Edge intersect test
    std::vector<DGUtils::Vec<3>> intersect_pts;
    intersecting_pts(s_c, nodeX_c, nodeY_c, nodeZ_c, intersect_pts);

    if(intersect_pts.size() == 3) {
      intersect_3pts(intersect_pts, sampleX_c, sampleY_c, sampleZ_c);
    } else if(intersect_pts.size() == 4) {
      intersect_4pts(intersect_pts, sampleX_c, sampleY_c, sampleZ_c);
    } else {
      for(int i = 0; i < LS_SAMPLE_NP; i++) {
        sampleX_c[i] = ref_r_ptr[i];
        sampleY_c[i] = ref_s_ptr[i];
        sampleZ_c[i] = ref_t_ptr[i];
        rst2xyz(sampleX_c[i], sampleY_c[i], sampleZ_c[i], nodeX_c, nodeY_c, nodeZ_c);
      }
    }

    int poly_ind = cell2polyMap.at(cell);
    PolyApprox3D pol = polys[poly_ind];
    DG_FP off_x, off_y, off_z;
    pol.get_offsets(off_x, off_y, off_z);
    for(int p = 0; p < LS_SAMPLE_NP; p++) {
      sampleX_c[p] -= off_x;
      sampleY_c[p] -= off_y;
      sampleZ_c[p] -= off_z;
    }

    // Simplified Newton sampling
    for(int p = 0; p < LS_SAMPLE_NP; p++) {
      bool converged = false;
      DG_FP start_x = sampleX_c[p];
      DG_FP start_y = sampleY_c[p];
      DG_FP start_z = sampleZ_c[p];
      if(!isnan(sampleX_c[p]))
        converged = simplified_newton(sampleX_c[p], sampleY_c[p], sampleZ_c[p], pol, tol);
      
      DG_FP dist_travelled = (start_x - sampleX_c[p]) * (start_x - sampleX_c[p]) 
                           + (start_y - sampleY_c[p]) * (start_y - sampleY_c[p])
                           + (start_z - sampleZ_c[p]) * (start_z - sampleZ_c[p]);

      // Check if point has converged
      // || !pt_in_tetra(sampleX_c[p], sampleY_c[p], sampleZ_c[p], nodeX_c, nodeY_c, nodeZ_c)
      // || dist_travelled > 1.5 * 1.5 * h * h
      if(!converged || !pt_in_tetra(sampleX_c[p] + off_x, sampleY_c[p] + off_y, sampleZ_c[p] + off_z, nodeX_c, nodeY_c, nodeZ_c)) {
        // Try randomly placing start and rerunning 10 times
        bool rerun_converged = false;
        bool rerun_in_bounds = false;
        int rerun_counter = 0;
        while((!rerun_converged || !rerun_in_bounds) && rerun_counter < 500) {
          // Random start point
          sampleX_c[p] = dis(gen);
          sampleY_c[p] = dis(gen);
          sampleZ_c[p] = dis(gen);
          rst2xyz(sampleX_c[p], sampleY_c[p], sampleZ_c[p], nodeX_c, nodeY_c, nodeZ_c);
          while(!pt_in_tetra(sampleX_c[p], sampleY_c[p], sampleZ_c[p], nodeX_c, nodeY_c, nodeZ_c)) {
            sampleX_c[p] = dis(gen);
            sampleY_c[p] = dis(gen);
            sampleZ_c[p] = dis(gen);
            rst2xyz(sampleX_c[p], sampleY_c[p], sampleZ_c[p], nodeX_c, nodeY_c, nodeZ_c);
          }
          sampleX_c[p] -= off_x;
          sampleY_c[p] -= off_y;
          sampleZ_c[p] -= off_z;
          // Rerun
          start_x = sampleX_c[p];
          start_y = sampleY_c[p];
          start_z = sampleZ_c[p];
          rerun_converged = simplified_newton(sampleX_c[p], sampleY_c[p], sampleZ_c[p], pol, tol);
          dist_travelled = (start_x - sampleX_c[p]) * (start_x - sampleX_c[p]) 
                         + (start_y - sampleY_c[p]) * (start_y - sampleY_c[p])
                         + (start_z - sampleZ_c[p]) * (start_z - sampleZ_c[p]);
          // rerun_in_bounds = dist_travelled <= 1.5 * 1.5 * h * h;
          rerun_in_bounds = pt_in_tetra(sampleX_c[p] + off_x, sampleY_c[p] + off_y, sampleZ_c[p] + off_z, nodeX_c, nodeY_c, nodeZ_c);
          rerun_counter++;
        }

        // || !rerun_in_bounds
        if(!rerun_converged || !rerun_in_bounds) {
          sampleX_c[p] = NAN;
          sampleY_c[p] = NAN;
          sampleZ_c[p] = NAN;
        }
      }
    }

    for(int p = 0; p < LS_SAMPLE_NP; p++) {
      sampleX_c[p] += off_x;
      sampleY_c[p] += off_y;
      sampleZ_c[p] += off_z;
    }
  }

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHost(mesh->z, OP_READ, z_ptr);

  releaseOP2PtrHost(s, OP_READ, s_ptr);
  releaseOP2PtrHost(mesh->nodeX, OP_READ, nodeX_ptr);
  releaseOP2PtrHost(mesh->nodeY, OP_READ, nodeY_ptr);
  releaseOP2PtrHost(mesh->nodeZ, OP_READ, nodeZ_ptr);
  releaseOP2PtrHost(sampleX, OP_WRITE, sampleX_ptr);
  releaseOP2PtrHost(sampleY, OP_WRITE, sampleY_ptr);
  releaseOP2PtrHost(sampleZ, OP_WRITE, sampleZ_ptr);
}
