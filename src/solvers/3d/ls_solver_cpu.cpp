#include "solvers/3d/ls_solver.h"

#include <set>

#include "dg_op2_blas.h"
#include "dg_utils.h"
#include "utils.h"
#include "ls_utils/3d/poly_approx.h"
#include "dg_dat_pool.h"

extern DGDatPool *dg_dat_pool;

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
/*
void set_sample_start_coords(DG_FP *r, DG_FP *s, DG_FP *t) {
  arma::vec x_, y_, z_, r_, s_, t_;
  DG3DUtils::setRefXYZ(2, x_, y_, z_);
  DG3DUtils::xyz2rst(x_, y_, z_, r_, s_, t_);
  for(int i = 0; i < LS_SAMPLE_NP; i++) {
    r[i] = r_[i];
    s[i] = s_[i];
    t[i] = t_[i];
  }
}
*/

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
  r[1] = (1.0 / 3.0) * (node0[0] - r_c) + r_c;
  s[1] = (1.0 / 3.0) * (node0[1] - s_c) + s_c;
  t[1] = (1.0 / 3.0) * (node0[2] - t_c) + t_c;
  r[2] = (2.0 / 3.0) * (node0[0] - r_c) + r_c;
  s[2] = (2.0 / 3.0) * (node0[1] - s_c) + s_c;
  t[2] = (2.0 / 3.0) * (node0[2] - t_c) + t_c;

  // Line from centre to node1
  r[3] = (1.0 / 3.0) * (node1[0] - r_c) + r_c;
  s[3] = (1.0 / 3.0) * (node1[1] - s_c) + s_c;
  t[3] = (1.0 / 3.0) * (node1[2] - t_c) + t_c;
  r[4] = (2.0 / 3.0) * (node1[0] - r_c) + r_c;
  s[4] = (2.0 / 3.0) * (node1[1] - s_c) + s_c;
  t[4] = (2.0 / 3.0) * (node1[2] - t_c) + t_c;

  // Line from centre to node2
  r[5] = (1.0 / 3.0) * (node2[0] - r_c) + r_c;
  s[5] = (1.0 / 3.0) * (node2[1] - s_c) + s_c;
  t[5] = (1.0 / 3.0) * (node2[2] - t_c) + t_c;
  r[6] = (2.0 / 3.0) * (node2[0] - r_c) + r_c;
  s[6] = (2.0 / 3.0) * (node2[1] - s_c) + s_c;
  t[6] = (2.0 / 3.0) * (node2[2] - t_c) + t_c;

  // Line from centre to node3
  r[7] = (1.0 / 3.0) * (node3[0] - r_c) + r_c;
  s[7] = (1.0 / 3.0) * (node3[1] - s_c) + s_c;
  t[7] = (1.0 / 3.0) * (node3[2] - t_c) + t_c;
  r[8] = (2.0 / 3.0) * (node3[0] - r_c) + r_c;
  s[8] = (2.0 / 3.0) * (node3[1] - s_c) + s_c;
  t[8] = (2.0 / 3.0) * (node3[2] - t_c) + t_c;

  r[9] = node0[0];
  s[9] = node0[1];
  t[9] = node0[2];
  r[10] = node1[0];
  s[10] = node1[1];
  t[10] = node1[2];
  r[11] = node2[0];
  s[11] = node2[1];
  t[11] = node2[2];
  r[12] = node3[0];
  s[12] = node3[1];
  t[12] = node3[2];
}

void intersect_3pts(const DG_FP *s, const DG_FP *nodeX, const DG_FP *nodeY,
                    const DG_FP *nodeZ, DG_FP *sampleX, DG_FP *sampleY,
                    DG_FP *sampleZ) {
  for(int i = 0; i < LS_SAMPLE_NP; i++) {
    sampleX[i] = NAN;
    sampleY[i] = NAN;
    sampleZ[i] = NAN;
  }
  // Get 3 intersecting points
  DG_FP inter_pts[3][3];
  int pt = 0;

  const DG_FP node0_s = s[0];
  const DG_FP node1_s = s[3];
  const DG_FP node2_s = s[9];
  const DG_FP node3_s = s[19];
  // Node0 -> Node1
  if(node0_s > 0.0 != node1_s > 0.0) {
    const DG_FP s_factor = node0_s / (node0_s - node1_s);
    inter_pts[pt][0] = nodeX[0] - s_factor * (nodeX[0] - nodeX[1]);
    inter_pts[pt][1] = nodeY[0] - s_factor * (nodeY[0] - nodeY[1]);
    inter_pts[pt][2] = nodeZ[0] - s_factor * (nodeZ[0] - nodeZ[1]);
    pt++;
  }
  // Node0 -> Node2
  if(node0_s > 0.0 != node2_s > 0.0) {
    const DG_FP s_factor = node0_s / (node0_s - node2_s);
    inter_pts[pt][0] = nodeX[0] - s_factor * (nodeX[0] - nodeX[2]);
    inter_pts[pt][1] = nodeY[0] - s_factor * (nodeY[0] - nodeY[2]);
    inter_pts[pt][2] = nodeZ[0] - s_factor * (nodeZ[0] - nodeZ[2]);
    pt++;
  }
  // Node0 -> Node3
  if(node0_s > 0.0 != node3_s > 0.0) {
    const DG_FP s_factor = node0_s / (node0_s - node3_s);
    inter_pts[pt][0] = nodeX[0] - s_factor * (nodeX[0] - nodeX[3]);
    inter_pts[pt][1] = nodeY[0] - s_factor * (nodeY[0] - nodeY[3]);
    inter_pts[pt][2] = nodeZ[0] - s_factor * (nodeZ[0] - nodeZ[3]);
    pt++;
  }
  // Node1 -> Node2
  if(node1_s > 0.0 != node2_s > 0.0) {
    const DG_FP s_factor = node1_s / (node1_s - node2_s);
    inter_pts[pt][0] = nodeX[1] - s_factor * (nodeX[1] - nodeX[2]);
    inter_pts[pt][1] = nodeY[1] - s_factor * (nodeY[1] - nodeY[2]);
    inter_pts[pt][2] = nodeZ[1] - s_factor * (nodeZ[1] - nodeZ[2]);
    pt++;
  }
  // Node1 -> Node3
  if(node1_s > 0.0 != node3_s > 0.0) {
    const DG_FP s_factor = node1_s / (node1_s - node3_s);
    inter_pts[pt][0] = nodeX[1] - s_factor * (nodeX[1] - nodeX[3]);
    inter_pts[pt][1] = nodeY[1] - s_factor * (nodeY[1] - nodeY[3]);
    inter_pts[pt][2] = nodeZ[1] - s_factor * (nodeZ[1] - nodeZ[3]);
    pt++;
  }
  // Node2 -> Node3
  if(node2_s > 0.0 != node3_s > 0.0) {
    const DG_FP s_factor = node2_s / (node2_s - node3_s);
    inter_pts[pt][0] = nodeX[2] - s_factor * (nodeX[2] - nodeX[3]);
    inter_pts[pt][1] = nodeY[2] - s_factor * (nodeY[2] - nodeY[3]);
    inter_pts[pt][2] = nodeZ[2] - s_factor * (nodeZ[2] - nodeZ[3]);
    pt++;
  }

  // Get intersecting plane eqn
  DG_FP plane_coeff[4];
  DG_FP vec0[3], vec1[3];
  vec0[0] = inter_pts[1][0] - inter_pts[0][0];
  vec0[1] = inter_pts[1][1] - inter_pts[0][1];
  vec0[2] = inter_pts[1][2] - inter_pts[0][2];
  vec1[0] = inter_pts[2][0] - inter_pts[0][0];
  vec1[1] = inter_pts[2][1] - inter_pts[0][1];
  vec1[2] = inter_pts[2][2] - inter_pts[0][2];
  plane_coeff[0] = vec0[1] * vec1[2] - vec0[2] * vec1[1];
  plane_coeff[1] = vec0[2] * vec1[0] - vec0[0] * vec1[2];
  plane_coeff[2] = vec0[0] * vec1[1] - vec0[1] * vec1[0];
  plane_coeff[3] = -(plane_coeff[0] * inter_pts[0][0] + plane_coeff[1] * inter_pts[0][1] + plane_coeff[2] * inter_pts[0][2]);

  // Calc sampling points
  DG_FP r_[] = {
    -1,
    -0.333333333333333,
    0.333333333333333,
    1,
    -1,
    -0.333333333333333,
    0.333333333333333,
    -1,
    -0.333333333333333,
    -1
  };
  DG_FP s_[] = {
    -1,
    -1,
    -1,
    -1,
    -0.333333333333333,
    -0.333333333333333,
    -0.333333333333333,
    0.333333333333333,
    0.333333333333333,
    1
  };

  DG_FP x_2D[10], y_2D[10];
  for(int i = 0; i < 10; i++) {
    x_2D[i] = 0.5*(-(r_[i]+s_[i])*inter_pts[0][0]+(1+r_[i])*inter_pts[1][0]+(1+s_[i])*inter_pts[2][0]);
    y_2D[i] = 0.5*(-(r_[i]+s_[i])*inter_pts[0][1]+(1+r_[i])*inter_pts[1][1]+(1+s_[i])*inter_pts[2][1]);
  }

  // Now project onto intersecting 3D plane
  for(int i = 0; i < 10; i++) {
    sampleX[i] = x_2D[i];
    sampleY[i] = y_2D[i];
    if(fabs(plane_coeff[2]) > 1e-8) {
      sampleZ[i] = (-plane_coeff[0] * x_2D[i] - plane_coeff[1] * y_2D[i] - plane_coeff[3]) / plane_coeff[2];
    } else {
      sampleX[i] = NAN;
      sampleY[i] = NAN;
      sampleZ[i] = NAN;
    }
  }
}

void intersect_4pts(const DG_FP *s, const DG_FP *nodeX, const DG_FP *nodeY,
                    const DG_FP *nodeZ, DG_FP *sampleX, DG_FP *sampleY,
                    DG_FP *sampleZ) {
  for(int i = 0; i < LS_SAMPLE_NP; i++) {
    sampleX[i] = NAN;
    sampleY[i] = NAN;
    sampleZ[i] = NAN;
  }
  // Get 3 intersecting points
  DG_FP inter_pts[4][3];
  int pt = 0;

  const DG_FP node0_s = s[0];
  const DG_FP node1_s = s[3];
  const DG_FP node2_s = s[9];
  const DG_FP node3_s = s[19];
  // Node0 -> Node1
  if(node0_s > 0.0 != node1_s > 0.0) {
    const DG_FP s_factor = node0_s / (node0_s - node1_s);
    inter_pts[pt][0] = nodeX[0] - s_factor * (nodeX[0] - nodeX[1]);
    inter_pts[pt][1] = nodeY[0] - s_factor * (nodeY[0] - nodeY[1]);
    inter_pts[pt][2] = nodeZ[0] - s_factor * (nodeZ[0] - nodeZ[1]);
    pt++;
  }
  // Node0 -> Node2
  if(node0_s > 0.0 != node2_s > 0.0) {
    const DG_FP s_factor = node0_s / (node0_s - node2_s);
    inter_pts[pt][0] = nodeX[0] - s_factor * (nodeX[0] - nodeX[2]);
    inter_pts[pt][1] = nodeY[0] - s_factor * (nodeY[0] - nodeY[2]);
    inter_pts[pt][2] = nodeZ[0] - s_factor * (nodeZ[0] - nodeZ[2]);
    pt++;
  }
  // Node0 -> Node3
  if(node0_s > 0.0 != node3_s > 0.0) {
    const DG_FP s_factor = node0_s / (node0_s - node3_s);
    inter_pts[pt][0] = nodeX[0] - s_factor * (nodeX[0] - nodeX[3]);
    inter_pts[pt][1] = nodeY[0] - s_factor * (nodeY[0] - nodeY[3]);
    inter_pts[pt][2] = nodeZ[0] - s_factor * (nodeZ[0] - nodeZ[3]);
    pt++;
  }
  // Node1 -> Node2
  if(node1_s > 0.0 != node2_s > 0.0) {
    const DG_FP s_factor = node1_s / (node1_s - node2_s);
    inter_pts[pt][0] = nodeX[1] - s_factor * (nodeX[1] - nodeX[2]);
    inter_pts[pt][1] = nodeY[1] - s_factor * (nodeY[1] - nodeY[2]);
    inter_pts[pt][2] = nodeZ[1] - s_factor * (nodeZ[1] - nodeZ[2]);
    pt++;
  }
  // Node1 -> Node3
  if(node1_s > 0.0 != node3_s > 0.0) {
    const DG_FP s_factor = node1_s / (node1_s - node3_s);
    inter_pts[pt][0] = nodeX[1] - s_factor * (nodeX[1] - nodeX[3]);
    inter_pts[pt][1] = nodeY[1] - s_factor * (nodeY[1] - nodeY[3]);
    inter_pts[pt][2] = nodeZ[1] - s_factor * (nodeZ[1] - nodeZ[3]);
    pt++;
  }
  // Node2 -> Node3
  if(node2_s > 0.0 != node3_s > 0.0) {
    const DG_FP s_factor = node2_s / (node2_s - node3_s);
    inter_pts[pt][0] = nodeX[2] - s_factor * (nodeX[2] - nodeX[3]);
    inter_pts[pt][1] = nodeY[2] - s_factor * (nodeY[2] - nodeY[3]);
    inter_pts[pt][2] = nodeZ[2] - s_factor * (nodeZ[2] - nodeZ[3]);
    pt++;
  }

  int edges[2][2];
  edges[0][0] = 0; edges[0][1] = 1; edges[1][0] = 2; edges[1][1] = 3;
  // Find closes point to point 0
  DG_FP min_dist = (inter_pts[0][0] - inter_pts[1][0]) * (inter_pts[0][0] - inter_pts[1][0])
                   + (inter_pts[0][1] - inter_pts[1][1]) * (inter_pts[0][1] - inter_pts[1][1])
                   + (inter_pts[0][2] - inter_pts[1][2]) * (inter_pts[0][2] - inter_pts[1][2]);
  DG_FP curr_dist = (inter_pts[0][0] - inter_pts[2][0]) * (inter_pts[0][0] - inter_pts[2][0])
                    + (inter_pts[0][1] - inter_pts[2][1]) * (inter_pts[0][1] - inter_pts[2][1])
                    + (inter_pts[0][2] - inter_pts[2][2]) * (inter_pts[0][2] - inter_pts[2][2]);
  if(curr_dist < min_dist) {
    min_dist = curr_dist;
    edges[0][1] = 2; edges[1][0] = 1; edges[1][1] = 3;
  }
  curr_dist = (inter_pts[0][0] - inter_pts[3][0]) * (inter_pts[0][0] - inter_pts[3][0])
             + (inter_pts[0][1] - inter_pts[3][1]) * (inter_pts[0][1] - inter_pts[3][1])
             + (inter_pts[0][2] - inter_pts[3][2]) * (inter_pts[0][2] - inter_pts[3][2]);
  if(curr_dist < min_dist) {
    min_dist = curr_dist;
    edges[0][1] = 3; edges[1][0] = 1; edges[1][1] = 2;
  }

  // Check both edges are going in the same direction
  DG_FP dist0 = (inter_pts[edges[0][0]][0] - inter_pts[edges[1][0]][0]) * (inter_pts[edges[0][0]][0] - inter_pts[edges[1][0]][0])
                + (inter_pts[edges[0][0]][1] - inter_pts[edges[1][0]][1]) * (inter_pts[edges[0][0]][1] - inter_pts[edges[1][0]][1])
                + (inter_pts[edges[0][0]][2] - inter_pts[edges[1][0]][2]) * (inter_pts[edges[0][0]][2] - inter_pts[edges[1][0]][2]);
  DG_FP dist1 = (inter_pts[edges[0][0]][0] - inter_pts[edges[1][1]][0]) * (inter_pts[edges[0][0]][0] - inter_pts[edges[1][1]][0])
                + (inter_pts[edges[0][0]][1] - inter_pts[edges[1][1]][1]) * (inter_pts[edges[0][0]][1] - inter_pts[edges[1][1]][1])
                + (inter_pts[edges[0][0]][2] - inter_pts[edges[1][1]][2]) * (inter_pts[edges[0][0]][2] - inter_pts[edges[1][1]][2]);
  if(dist1 < dist0) {
    int tmp = edges[1][0];
    edges[1][0] = edges[1][1];
    edges[1][1] = tmp;
  }

  // Have intersection quad, now place sample points
  for(int i = 0; i < 3; i++) {
    DG_FP f = (1.0 / 2.0) * i;
    DG_FP pt_0[3];
    pt_0[0] = f * (inter_pts[edges[0][1]][0] - inter_pts[edges[0][0]][0]) + inter_pts[edges[0][0]][0];
    pt_0[1] = f * (inter_pts[edges[0][1]][1] - inter_pts[edges[0][0]][1]) + inter_pts[edges[0][0]][1];
    pt_0[2] = f * (inter_pts[edges[0][1]][2] - inter_pts[edges[0][0]][2]) + inter_pts[edges[0][0]][2];

    DG_FP pt_1[3];
    pt_1[0] = f * (inter_pts[edges[1][1]][0] - inter_pts[edges[1][0]][0]) + inter_pts[edges[1][0]][0];
    pt_1[1] = f * (inter_pts[edges[1][1]][1] - inter_pts[edges[1][0]][1]) + inter_pts[edges[1][0]][1];
    pt_1[2] = f * (inter_pts[edges[1][1]][2] - inter_pts[edges[1][0]][2]) + inter_pts[edges[1][0]][2];

    for(int j = 0; j < 3; j++) {
      DG_FP fj = (1.0 / 2.0) * j;
      sampleX[i * 3 + j] = fj * (pt_1[0] - pt_0[0]) + pt_0[0];
      sampleY[i * 3 + j] = fj * (pt_1[1] - pt_0[1]) + pt_0[1];
      sampleZ[i * 3 + j] = fj * (pt_1[2] - pt_0[2]) + pt_0[2];
    }
  }
}


void LevelSetSolver3D::sampleInterface(op_dat sampleX, op_dat sampleY, op_dat sampleZ) {
  DG_FP ref_r[LS_SAMPLE_NP], ref_s[LS_SAMPLE_NP], ref_t[LS_SAMPLE_NP];
  set_sample_start_coords(ref_r, ref_s, ref_t);
  DGTempDat tmp_s_modal = dg_dat_pool->requestTempDatCells(DG_NP);
  op2_gemv(mesh, false, 1.0, DGConstants::INV_V, s, 0.0, tmp_s_modal.dat);

  const DG_FP *ref_r_ptr = ref_r;
  const DG_FP *ref_s_ptr = ref_s;
  const DG_FP *ref_t_ptr = ref_t;
  const DG_FP *s_ptr = getOP2PtrHost(s, OP_READ);
  const DG_FP *s_modal_ptr = getOP2PtrHost(tmp_s_modal.dat, OP_READ);
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
    const DG_FP *s_c = s_ptr + cell * DG_NP;
    const DG_FP *s_modal_c = s_modal_ptr + cell * DG_NP;
    const DG_FP *nodeX_c = nodeX_ptr + cell * 4;
    const DG_FP *nodeY_c = nodeY_ptr + cell * 4;
    const DG_FP *nodeZ_c = nodeZ_ptr + cell * 4;
    const DG_FP *x_c = x_ptr + cell * DG_NP;
    const DG_FP *y_c = y_ptr + cell * DG_NP;
    const DG_FP *z_c = z_ptr + cell * DG_NP;
    DG_FP *sampleX_c = sampleX_ptr + cell * LS_SAMPLE_NP;
    DG_FP *sampleY_c = sampleY_ptr + cell * LS_SAMPLE_NP;
    DG_FP *sampleZ_c = sampleZ_ptr + cell * LS_SAMPLE_NP;

    bool positive0 = s_c[0] > 0.0;
    bool interface = false;
    for(int i = 1; i < DG_NP; i++) {
      if(positive0 != s_c[i] > 0.0)
        interface = true;
    }
    if(!interface) {
      for(int i = 0; i < LS_SAMPLE_NP; i++) {
        sampleX_c[i] = NAN;
        sampleY_c[i] = NAN;
        sampleZ_c[i] = NAN;
      }
      continue;
    }

    // Edge intersect test
    int edge_count = 0;
    const DG_FP node0_s = s_c[0];
    const DG_FP node1_s = s_c[3];
    const DG_FP node2_s = s_c[9];
    const DG_FP node3_s = s_c[19];
    // Node0 -> Node1
    if(node0_s > 0.0 != node1_s > 0.0)
      edge_count++;
    // Node0 -> Node2
    if(node0_s > 0.0 != node2_s > 0.0)
      edge_count++;
    // Node0 -> Node3
    if(node0_s > 0.0 != node3_s > 0.0)
      edge_count++;
    // Node1 -> Node2
    if(node1_s > 0.0 != node2_s > 0.0)
      edge_count++;
    // Node1 -> Node3
    if(node1_s > 0.0 != node3_s > 0.0)
      edge_count++;
    // Node2 -> Node3
    if(node2_s > 0.0 != node3_s > 0.0)
      edge_count++;
    // End

    if(edge_count == 3) {
      intersect_3pts(s_c, nodeX_c, nodeY_c, nodeZ_c, sampleX_c, sampleY_c, sampleZ_c);
    } else if(edge_count == 4) {
      intersect_4pts(s_c, nodeX_c, nodeY_c, nodeZ_c, sampleX_c, sampleY_c, sampleZ_c);
    } else {
      // Revert back to simplified Newton sampling
      bool smpPt = false;
      for(int i = 0; i < LS_SAMPLE_NP; i++) {
        sampleX_c[i] = ref_r_ptr[i];
        sampleY_c[i] = ref_s_ptr[i];
        sampleZ_c[i] = ref_t_ptr[i];
      }

      for(int p = 0; p < LS_SAMPLE_NP; p++) {
        bool converged = false;
        for(int step = 0; step < 10; step++) {
          DG_FP surf = DGUtils::val_at_pt_3d(sampleX_c[p], sampleY_c[p], sampleZ_c[p], s_modal_c, DG_ORDER);
          DG_FP dsdx, dsdy, dsdz;
          DGUtils::grad_at_pt_3d(sampleX_c[p], sampleY_c[p], sampleZ_c[p], s_modal_c, DG_ORDER, dsdx, dsdy, dsdz);

          DG_FP sqrnorm = dsdx * dsdx + dsdy * dsdy + dsdz * dsdz;
          if(sqrnorm > 0.0) {
            dsdx *= surf / sqrnorm;
            dsdy *= surf / sqrnorm;
            dsdz *= surf / sqrnorm;
          }

          sampleX_c[p] -= dsdx;
          sampleY_c[p] -= dsdy;
          sampleZ_c[p] -= dsdz;

          // Check convergence
          if(dsdx * dsdx + dsdy * dsdy + dsdz * dsdz < 1e-5) {
            converged = true;
            break;
          }
        }

        // Convert back to xyz coords
        if(converged) {
          rst2xyz(sampleX_c[p], sampleY_c[p], sampleZ_c[p], nodeX_c, nodeY_c, nodeZ_c);

          // Check if point is in cell
          if(!pt_in_tetra(sampleX_c[p], sampleY_c[p], sampleZ_c[p], nodeX_c, nodeY_c, nodeZ_c)) {
            sampleX_c[p] = NAN;
            sampleY_c[p] = NAN;
            sampleZ_c[p] = NAN;
          } else {
            smpPt = true;
          }
        } else {
          sampleX_c[p] = NAN;
          sampleY_c[p] = NAN;
          sampleZ_c[p] = NAN;
        }
      }

      if(!smpPt) {
        int id = 0;
        DG_FP val = fabs(s_c[0]);
        for(int i = 1; i < DG_NP; i++) {
          if(fabs(s_c[i]) < val) {
            val = fabs(s_c[i]);
            id = i;
          }
        }

        sampleX_c[0] = x_c[id];
        sampleY_c[0] = y_c[id];
        sampleZ_c[0] = z_c[id];
      }
    }
  }

  releaseOP2PtrHost(mesh->x, OP_READ, x_ptr);
  releaseOP2PtrHost(mesh->y, OP_READ, y_ptr);
  releaseOP2PtrHost(mesh->z, OP_READ, z_ptr);

  releaseOP2PtrHost(s, OP_READ, s_ptr);
  releaseOP2PtrHost(tmp_s_modal.dat, OP_READ, s_modal_ptr);
  releaseOP2PtrHost(mesh->nodeX, OP_READ, nodeX_ptr);
  releaseOP2PtrHost(mesh->nodeY, OP_READ, nodeY_ptr);
  releaseOP2PtrHost(mesh->nodeZ, OP_READ, nodeZ_ptr);
  releaseOP2PtrHost(sampleX, OP_WRITE, sampleX_ptr);
  releaseOP2PtrHost(sampleY, OP_WRITE, sampleY_ptr);
  releaseOP2PtrHost(sampleZ, OP_WRITE, sampleZ_ptr);

  dg_dat_pool->releaseTempDatCells(tmp_s_modal);
}
