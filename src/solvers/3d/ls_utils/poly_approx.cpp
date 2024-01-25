#include "ls_utils/3d/poly_approx.h"

#include <random>

#include "timing.h"
#include "dg_constants/dg_constants_3d.h"
#include "dg_utils.h"

extern Timing *timer;
extern DGConstants *constants;

using namespace std;

PolyApprox3D::PolyApprox3D(const int cell_ind, set<int> stencil,
                       const DG_FP *x_ptr, const DG_FP *y_ptr,
                       const DG_FP *z_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr) {
  get_offset(cell_ind, x_ptr, y_ptr, z_ptr);

  vector<DG_FP> x_vec, y_vec, z_vec, s_vec;
  stencil_data(cell_ind, stencil, x_ptr, y_ptr, z_ptr, s_ptr, modal_ptr, x_vec, y_vec, z_vec, s_vec);

  if(N == 2) {
    set_2nd_order_coeff(x_vec, y_vec, z_vec, s_vec);
  } else if(N == 3) {
    set_3rd_order_coeff(x_vec, y_vec, z_vec, s_vec);
  } else if(N == 4) {
    set_4th_order_coeff(x_vec, y_vec, z_vec, s_vec);
  }
}

PolyApprox3D::PolyApprox3D(std::vector<DG_FP> &c, DG_FP off_x, DG_FP off_y,
                       DG_FP off_z) {
  offset_x = off_x;
  offset_y = off_y;
  offset_z = off_z;
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(c[i]);
  }
}

void PolyApprox3D::get_offset(const int ind, const DG_FP *x_ptr, const DG_FP *y_ptr,
                            const DG_FP *z_ptr) {
  offset_x = x_ptr[ind * DG_NP];
  offset_y = y_ptr[ind * DG_NP];
  offset_z = z_ptr[ind * DG_NP];
  // offset_x = 0.0;
  // offset_y = 0.0;
  // offset_z = 0.0;
}

struct Coord {
  DG_FP x;
  DG_FP y;
  DG_FP z;
};

struct Point {
  Coord coord;
  DG_FP val;
  int count;
};

struct cmpCoords {
    bool operator()(const Coord& a, const Coord& b) const {
        bool xCmp = abs(a.x - b.x) < 1e-8;
        bool yCmp = abs(a.y - b.y) < 1e-8;
        bool zCmp = abs(a.z - b.z) < 1e-8;
        if(xCmp && yCmp && zCmp) {
          return false;
        } else if(xCmp && yCmp) {
          return a.z < b.z;
        } else if(xCmp) {
          return a.y < b.y;
        } else {
          return a.x < b.x;
        }
    }
};

bool in_tetra(const DG_FP ptX, const DG_FP ptY, const DG_FP ptZ,
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

void pol_rst2xyz(DG_FP &sampleX, DG_FP &sampleY, DG_FP &sampleZ,
             const DG_FP *nodeX, const DG_FP *nodeY,
             const DG_FP *nodeZ) {
  DG_FP r_ = sampleX;
  DG_FP s_ = sampleY;
  DG_FP t_ = sampleZ;

  sampleX = 0.5 * (-(1.0 + r_ + s_ + t_) * nodeX[0] + (1.0 + r_) * nodeX[1] + (1.0 + s_) * nodeX[2] + (1.0 + t_) * nodeX[3]);
  sampleY = 0.5 * (-(1.0 + r_ + s_ + t_) * nodeY[0] + (1.0 + r_) * nodeY[1] + (1.0 + s_) * nodeY[2] + (1.0 + t_) * nodeY[3]);
  sampleZ = 0.5 * (-(1.0 + r_ + s_ + t_) * nodeZ[0] + (1.0 + r_) * nodeZ[1] + (1.0 + s_) * nodeZ[2] + (1.0 + t_) * nodeZ[3]);
}

bool pol_simplified_newton(DG_FP &pt_r, DG_FP &pt_s, DG_FP &pt_t, const DG_FP *modal,
                    const DG_FP tol, const DG_FP *tetra_r, const DG_FP *tetra_s,
                    const DG_FP *tetra_t) {
  bool converged = false;
  for(int step = 0; step < 20; step++) {
    DG_FP surf = DGUtils::val_at_pt_3d(pt_r, pt_s, pt_t, modal, DG_ORDER);
    DG_FP dsdr, dsds, dsdt;
    DGUtils::grad_at_pt_3d(pt_r, pt_s, pt_t, modal, DG_ORDER, dsdr, dsds, dsdt);

    DG_FP sqrnorm = dsdr * dsdr + dsds * dsds + dsdt * dsdt;
    if(sqrnorm > 0.0) {
      dsdr *= surf / sqrnorm;
      dsds *= surf / sqrnorm;
      dsdt *= surf / sqrnorm;
    }

    pt_r -= dsdr;
    pt_s -= dsds;
    pt_t -= dsdt;

    if(!in_tetra(pt_r, pt_s, pt_t, tetra_r, tetra_s, tetra_t)) {
      break;
    }

    // Check convergence
    if(dsdr * dsdr + dsds * dsds + dsdt * dsdt < tol) {
      converged = true;
      break;
    }
  }
  return converged;
}

void PolyApprox3D::stencil_data(const int cell_ind, const set<int> &stencil,
                              const DG_FP *x_ptr, const DG_FP *y_ptr,
                              const DG_FP *z_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr,
                              vector<DG_FP> &x, vector<DG_FP> &y,
                              vector<DG_FP> &z, vector<DG_FP> &s) {
  // Setup random number generator for later
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1.0, 1.0);

  map<Coord, Point, cmpCoords> pointMap;

  for(int n = 0; n < DG_NP; n++) {
    int ind = cell_ind * DG_NP + n;
    Coord coord;
    coord.x = x_ptr[ind] - offset_x;
    coord.y = y_ptr[ind] - offset_y;
    coord.z = z_ptr[ind] - offset_z;
    Point point;
    auto res = pointMap.insert(make_pair(coord, point));
    res.first->second.coord = coord;
    res.first->second.val   = s_ptr[ind];
    res.first->second.count = 1;
  }

/*
  std::vector<bool> consider;
  for(const auto &sten : stencil) {
    if(sten == cell_ind) {
      consider.push_back(false);
      continue;
    }

    bool share_pt = false;
    for(int n = 0; n < DG_NP; n++) {
      int ind = sten * DG_NP + n;
      Coord coord;
      coord.x = x_ptr[ind] - offset_x;
      coord.y = y_ptr[ind] - offset_y;
      coord.z = z_ptr[ind] - offset_z;
      Point point;
      auto res = pointMap.find(coord);
      if(res != pointMap.end()) {
        share_pt = true;
        break;
      }
    }
    consider.push_back(share_pt);
  }

  int i = 0;
  for(const auto &sten : stencil) {
    if(!consider[i]) {
      i++;
      continue;
    }

    i++;
    for(int n = 0; n < DG_NP; n++) {
      int ind = sten * DG_NP + n;

      Coord coord;
      coord.x = x_ptr[ind] - offset_x;
      coord.y = y_ptr[ind] - offset_y;
      coord.z = z_ptr[ind] - offset_z;
      Point point;
      auto res = pointMap.insert(make_pair(coord, point));

      if(res.second) {
        // Point was inserted
        res.first->second.coord = coord;
        res.first->second.val   = s_ptr[ind];
        res.first->second.count = 1;
      } else {
        // Point already exists
        res.first->second.val += s_ptr[ind];
        res.first->second.count++;
      }
    }
  }
*/

  for(const auto &sten : stencil) {
    if(sten == cell_ind) continue;
    for(int n = 0; n < DG_NP; n++) {
      int ind = sten * DG_NP + n;

      Coord coord;
      coord.x = x_ptr[ind] - offset_x;
      coord.y = y_ptr[ind] - offset_y;
      coord.z = z_ptr[ind] - offset_z;
      Point point;
      auto res = pointMap.find(coord);
      if(res == pointMap.end()) continue;

      res->second.val += s_ptr[ind];
      res->second.count++;
    }
  }

/*
  for(const auto &sten : stencil) {
    if(sten != cell_ind) continue;
    for(int n = 0; n < DG_NP; n++) {
      int ind = sten * DG_NP + n;

      Coord coord;
      coord.x = x_ptr[ind] - offset_x;
      coord.y = y_ptr[ind] - offset_y;
      coord.z = z_ptr[ind] - offset_z;
      Point point;
      auto res = pointMap.insert(make_pair(coord, point));

      if(res.second) {
        // Point was inserted
        res.first->second.coord = coord;
        res.first->second.val   = s_ptr[ind];
        res.first->second.count = 1;
      } else {
        // if(sten == cell_ind) {
        //   res.first->second.val = s_ptr[ind];
        // }
        // Point already exists
        res.first->second.val += s_ptr[ind];
        res.first->second.count++;
      }
    }
  }
*/
  const DG_FP *tmp_r_ptr = constants->get_mat_ptr(DGConstants::R) + (DG_ORDER - 1) * DG_NP;
  const DG_FP *tmp_s_ptr = constants->get_mat_ptr(DGConstants::S) + (DG_ORDER - 1) * DG_NP;
  const DG_FP *tmp_t_ptr = constants->get_mat_ptr(DGConstants::T) + (DG_ORDER - 1) * DG_NP;

  const DG_FP tetra_r[4] = {tmp_r_ptr[0], tmp_r_ptr[3], tmp_r_ptr[9], tmp_r_ptr[19]};
  const DG_FP tetra_s[4] = {tmp_s_ptr[0], tmp_s_ptr[3], tmp_s_ptr[9], tmp_s_ptr[19]};
  const DG_FP tetra_t[4] = {tmp_t_ptr[0], tmp_t_ptr[3], tmp_t_ptr[9], tmp_t_ptr[19]};

  for(const auto &sten : stencil) {
    bool contain_interface = false;
    bool s_pos = s_ptr[sten * DG_NP] > 0.0;
    for(int i = 0; i < DG_NP; i++) {
      if(s_ptr[sten * DG_NP + i] > 0.0 != s_pos) {
        contain_interface = true;
        break;
      }
    }
    if(!contain_interface && sten == cell_ind) continue;

    DG_FP ref_r_ptr[DG_NP], ref_s_ptr[DG_NP], ref_t_ptr[DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      ref_r_ptr[i] = tmp_r_ptr[i] * 0.75;
      ref_s_ptr[i] = tmp_s_ptr[i] * 0.75;
      ref_t_ptr[i] = tmp_t_ptr[i] * 0.75;
    }

    const DG_FP X0 = x_ptr[sten * DG_NP + 0]; const DG_FP Y0 = y_ptr[sten * DG_NP + 0]; const DG_FP Z0 = z_ptr[sten * DG_NP + 0];
    const DG_FP X1 = x_ptr[sten * DG_NP + 3]; const DG_FP Y1 = y_ptr[sten * DG_NP + 3]; const DG_FP Z1 = z_ptr[sten * DG_NP + 3];
    const DG_FP X2 = x_ptr[sten * DG_NP + 9]; const DG_FP Y2 = y_ptr[sten * DG_NP + 9]; const DG_FP Z2 = z_ptr[sten * DG_NP + 9];
    const DG_FP X3 = x_ptr[sten * DG_NP + 19]; const DG_FP Y3 = y_ptr[sten * DG_NP + 19]; const DG_FP Z3 = z_ptr[sten * DG_NP + 19];

    #pragma omp parallel for
    for(int p = 0; p < DG_NP; p++) {
      // bool converged = false;
      bool converged = pol_simplified_newton(ref_r_ptr[p], ref_s_ptr[p], ref_t_ptr[p], &modal_ptr[sten * DG_NP], 1e-8, tetra_r, tetra_s, tetra_t);

      // if(!converged && sten == cell_ind) {
      //   int counter = 0;
      //   while(!converged && counter < 10) {
      //     ref_r_ptr[p] = dis(gen);
      //     ref_s_ptr[p] = dis(gen);
      //     ref_t_ptr[p] = dis(gen);
      //     while(!in_tetra(ref_r_ptr[p], ref_s_ptr[p], ref_t_ptr[p], tetra_r, tetra_s, tetra_t)) {
      //       ref_r_ptr[p] = dis(gen);
      //       ref_s_ptr[p] = dis(gen);
      //       ref_t_ptr[p] = dis(gen);
      //     }
      //     converged = pol_simplified_newton(ref_r_ptr[p], ref_s_ptr[p], ref_t_ptr[p], &modal_ptr[sten * DG_NP], 1e-8, tetra_r, tetra_s, tetra_t);
      //     counter++;
      //   }
      // }

      // Check if point has converged
      // && in_tetra(ref_r_ptr[p], ref_s_ptr[p], ref_t_ptr[p], tetra_r, tetra_s, tetra_t)
      if(converged) {
        const DG_FP x = 0.5 * (-(1.0 + ref_r_ptr[p] + ref_s_ptr[p] + ref_t_ptr[p]) * X0 + (1.0 + ref_r_ptr[p]) * X1 + (1.0 + ref_s_ptr[p]) * X2 + (1.0 + ref_t_ptr[p]) * X3);
        const DG_FP y = 0.5 * (-(1.0 + ref_r_ptr[p] + ref_s_ptr[p] + ref_t_ptr[p]) * Y0 + (1.0 + ref_r_ptr[p]) * Y1 + (1.0 + ref_s_ptr[p]) * Y2 + (1.0 + ref_t_ptr[p]) * Y3);
        const DG_FP z = 0.5 * (-(1.0 + ref_r_ptr[p] + ref_s_ptr[p] + ref_t_ptr[p]) * Z0 + (1.0 + ref_r_ptr[p]) * Z1 + (1.0 + ref_s_ptr[p]) * Z2 + (1.0 + ref_t_ptr[p]) * Z3);
        const DG_FP tmp_val = DGUtils::val_at_pt_3d(ref_r_ptr[p], ref_s_ptr[p], ref_t_ptr[p], &modal_ptr[sten * DG_NP], DG_ORDER);
        #pragma omp critical
        {
          Coord coord;
          coord.x = x - offset_x;
          coord.y = y - offset_y;
          coord.z = z - offset_z;
          Point point;
          auto res = pointMap.insert(make_pair(coord, point));
          if(res.second) {
            // Point was inserted
            res.first->second.coord = coord;
            res.first->second.val   = tmp_val;
            res.first->second.count = 1;
          } else {
            // Point already exists
            res.first->second.val += tmp_val;
            res.first->second.count++;
          }
        }
      }
    }
  }

  const DG_FP X0 = x_ptr[cell_ind * DG_NP + 0]; const DG_FP Y0 = y_ptr[cell_ind * DG_NP + 0]; const DG_FP Z0 = z_ptr[cell_ind * DG_NP + 0];
  const DG_FP X1 = x_ptr[cell_ind * DG_NP + 3]; const DG_FP Y1 = y_ptr[cell_ind * DG_NP + 3]; const DG_FP Z1 = z_ptr[cell_ind * DG_NP + 3];
  const DG_FP X2 = x_ptr[cell_ind * DG_NP + 9]; const DG_FP Y2 = y_ptr[cell_ind * DG_NP + 9]; const DG_FP Z2 = z_ptr[cell_ind * DG_NP + 9];
  const DG_FP X3 = x_ptr[cell_ind * DG_NP + 19]; const DG_FP Y3 = y_ptr[cell_ind * DG_NP + 19]; const DG_FP Z3 = z_ptr[cell_ind * DG_NP + 19];
  const DG_FP nodeX[] = {X0, X1, X2, X3};
  const DG_FP nodeY[] = {Y0, Y1, Y2, Y3};
  const DG_FP nodeZ[] = {Z0, Z1, Z2, Z3};

  for(int p = 0; p < DG_NP; p++) {
    DG_FP sampleX = dis(gen);
    DG_FP sampleY = dis(gen);
    DG_FP sampleZ = dis(gen);
    DG_FP surf = DGUtils::val_at_pt_3d(sampleX, sampleY, sampleZ, &modal_ptr[cell_ind * DG_NP], DG_ORDER);
    pol_rst2xyz(sampleX, sampleY, sampleZ, nodeX, nodeY, nodeZ);
    while(!in_tetra(sampleX, sampleY, sampleZ, nodeX, nodeY, nodeZ)) {
      sampleX = dis(gen);
      sampleY = dis(gen);
      sampleZ = dis(gen);
      surf = DGUtils::val_at_pt_3d(sampleX, sampleY, sampleZ, &modal_ptr[cell_ind * DG_NP], DG_ORDER);
      pol_rst2xyz(sampleX, sampleY, sampleZ, nodeX, nodeY, nodeZ);
    }

    Coord coord;
    coord.x = sampleX - offset_x;
    coord.y = sampleY - offset_y;
    coord.z = sampleZ - offset_z;
    Point point;
    auto res = pointMap.insert(make_pair(coord, point));
    if(res.second) {
      // Point was inserted
      res.first->second.coord = coord;
      res.first->second.val   = surf;
      res.first->second.count = 1;
    }
  }


  DG_FP ref_r_ptr[DG_NP], ref_s_ptr[DG_NP], ref_t_ptr[DG_NP];
  for(int i = 0; i < DG_NP; i++) {
    ref_r_ptr[i] = tmp_r_ptr[i] * 0.75;
    ref_s_ptr[i] = tmp_s_ptr[i] * 0.75;
    ref_t_ptr[i] = tmp_t_ptr[i] * 0.75;
  }
/*
  for(const auto &sten : stencil) {
    const DG_FP X0 = x_ptr[sten * DG_NP + 0]; const DG_FP Y0 = y_ptr[sten * DG_NP + 0]; const DG_FP Z0 = z_ptr[sten * DG_NP + 0];
    const DG_FP X1 = x_ptr[sten * DG_NP + 3]; const DG_FP Y1 = y_ptr[sten * DG_NP + 3]; const DG_FP Z1 = z_ptr[sten * DG_NP + 3];
    const DG_FP X2 = x_ptr[sten * DG_NP + 9]; const DG_FP Y2 = y_ptr[sten * DG_NP + 9]; const DG_FP Z2 = z_ptr[sten * DG_NP + 9];
    const DG_FP X3 = x_ptr[sten * DG_NP + 19]; const DG_FP Y3 = y_ptr[sten * DG_NP + 19]; const DG_FP Z3 = z_ptr[sten * DG_NP + 19];
    for(int n = 0; n < DG_NP; n++) {
      int ind = sten * DG_NP + n;

      const DG_FP r = ref_r_ptr[n];
      const DG_FP s = ref_s_ptr[n];
      const DG_FP t = ref_t_ptr[n];
      const DG_FP val = DGUtils::val_at_pt_3d(r, s, t, &modal_ptr[sten * DG_NP], DG_ORDER);

      const DG_FP x = 0.5 * (-(1.0 + r + s + t) * X0 + (1.0 + r) * X1 + (1.0 + s) * X2 + (1.0 + t) * X3);
      const DG_FP y = 0.5 * (-(1.0 + r + s + t) * Y0 + (1.0 + r) * Y1 + (1.0 + s) * Y2 + (1.0 + t) * Y3);
      const DG_FP z = 0.5 * (-(1.0 + r + s + t) * Z0 + (1.0 + r) * Z1 + (1.0 + s) * Z2 + (1.0 + t) * Z3);

      Coord coord;
      coord.x = x - offset_x;
      coord.y = y - offset_y;
      coord.z = z - offset_z;
      Point point;
      auto res = pointMap.insert(make_pair(coord, point));

      if(res.second) {
        // Point was inserted
        res.first->second.coord = coord;
        res.first->second.val   = val;
        res.first->second.count = 1;
      } else {
        // Point already exists
        res.first->second.val += val;
        res.first->second.count++;
      }
    }
  }
*/
  for(auto const &p : pointMap) {
    x.push_back(p.second.coord.x);
    y.push_back(p.second.coord.y);
    z.push_back(p.second.coord.z);
    s.push_back(p.second.val / (DG_FP)p.second.count);
  }
}

#define C_POLY_IND 0
#define X_POLY_IND 1
#define Y_POLY_IND 2
#define Z_POLY_IND 3
#define XY_POLY_IND 4
#define XZ_POLY_IND 5
#define YZ_POLY_IND 6
#define X2_POLY_IND 7
#define Y2_POLY_IND 8
#define Z2_POLY_IND 9
#define XYZ_POLY_IND 10
#define X2Y_POLY_IND 11
#define X2Z_POLY_IND 12
#define Y2X_POLY_IND 13
#define Y2Z_POLY_IND 14
#define Z2X_POLY_IND 15
#define Z2Y_POLY_IND 16
#define X3_POLY_IND 17
#define Y3_POLY_IND 18
#define Z3_POLY_IND 19
#define X2YZ_POLY_IND 20
#define XY2Z_POLY_IND 21
#define XYZ2_POLY_IND 22
#define X2Y2_POLY_IND 23
#define X2Z2_POLY_IND 24
#define Y2Z2_POLY_IND 25
#define X3Y_POLY_IND 26
#define X3Z_POLY_IND 27
#define Y3X_POLY_IND 28
#define Y3Z_POLY_IND 29
#define Z3X_POLY_IND 30
#define Z3Y_POLY_IND 31
#define X4_POLY_IND 32
#define Y4_POLY_IND 33
#define Z4_POLY_IND 34

void PolyApprox3D::set_2nd_order_coeff(const vector<DG_FP> &x, const vector<DG_FP> &y,
                                     const vector<DG_FP> &z, const vector<DG_FP> &s) {
  arma::mat A(x.size(), num_coeff());
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,C_POLY_IND)  = 1.0;
    A(i,X_POLY_IND)  = x[i];
    A(i,Y_POLY_IND)  = y[i];
    A(i,Z_POLY_IND)  = z[i];
    A(i,XY_POLY_IND) = x[i] * y[i];
    A(i,XZ_POLY_IND) = x[i] * z[i];
    A(i,YZ_POLY_IND) = y[i] * z[i];
    A(i,X2_POLY_IND) = x[i] * x[i];
    A(i,Y2_POLY_IND) = y[i] * y[i];
    A(i,Z2_POLY_IND) = z[i] * z[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(ans(i));
  }
}

void PolyApprox3D::set_3rd_order_coeff(const vector<DG_FP> &x, const vector<DG_FP> &y,
                                     const vector<DG_FP> &z, const vector<DG_FP> &s) {
  arma::mat A(x.size(), num_coeff());
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,C_POLY_IND)   = 1.0;
    A(i,X_POLY_IND)   = x[i];
    A(i,Y_POLY_IND)   = y[i];
    A(i,Z_POLY_IND)   = z[i];
    A(i,XY_POLY_IND)  = x[i] * y[i];
    A(i,XZ_POLY_IND)  = x[i] * z[i];
    A(i,YZ_POLY_IND)  = y[i] * z[i];
    A(i,X2_POLY_IND)  = x[i] * x[i];
    A(i,Y2_POLY_IND)  = y[i] * y[i];
    A(i,Z2_POLY_IND)  = z[i] * z[i];
    A(i,XYZ_POLY_IND) = x[i] * y[i] * z[i];
    A(i,X2Y_POLY_IND) = x[i] * x[i] * y[i];
    A(i,X2Z_POLY_IND) = x[i] * x[i] * z[i];
    A(i,Y2X_POLY_IND) = y[i] * y[i] * x[i];
    A(i,Y2Z_POLY_IND) = y[i] * y[i] * z[i];
    A(i,Z2X_POLY_IND) = z[i] * z[i] * x[i];
    A(i,Z2Y_POLY_IND) = z[i] * z[i] * y[i];
    A(i,X3_POLY_IND)  = x[i] * x[i] * x[i];
    A(i,Y3_POLY_IND)  = y[i] * y[i] * y[i];
    A(i,Z3_POLY_IND)  = z[i] * z[i] * z[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(ans(i));
  }
}

void PolyApprox3D::set_4th_order_coeff(const vector<DG_FP> &x, const vector<DG_FP> &y,
                                     const vector<DG_FP> &z, const vector<DG_FP> &s) {
  arma::mat A(x.size(), num_coeff());
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,C_POLY_IND)    = 1.0;
    A(i,X_POLY_IND)    = x[i];
    A(i,Y_POLY_IND)    = y[i];
    A(i,Z_POLY_IND)    = z[i];
    A(i,XY_POLY_IND)   = x[i] * y[i];
    A(i,XZ_POLY_IND)   = x[i] * z[i];
    A(i,YZ_POLY_IND)   = y[i] * z[i];
    A(i,X2_POLY_IND)   = x[i] * x[i];
    A(i,Y2_POLY_IND)   = y[i] * y[i];
    A(i,Z2_POLY_IND)   = z[i] * z[i];
    A(i,XYZ_POLY_IND)  = x[i] * y[i] * z[i];
    A(i,X2Y_POLY_IND)  = x[i] * x[i] * y[i];
    A(i,X2Z_POLY_IND)  = x[i] * x[i] * z[i];
    A(i,Y2X_POLY_IND)  = y[i] * y[i] * x[i];
    A(i,Y2Z_POLY_IND)  = y[i] * y[i] * z[i];
    A(i,Z2X_POLY_IND)  = z[i] * z[i] * x[i];
    A(i,Z2Y_POLY_IND)  = z[i] * z[i] * y[i];
    A(i,X3_POLY_IND)   = x[i] * x[i] * x[i];
    A(i,Y3_POLY_IND)   = y[i] * y[i] * y[i];
    A(i,Z3_POLY_IND)   = z[i] * z[i] * z[i];
    A(i,X2YZ_POLY_IND) = x[i] * x[i] * y[i] * z[i];
    A(i,XY2Z_POLY_IND) = x[i] * y[i] * y[i] * z[i];
    A(i,XYZ2_POLY_IND) = x[i] * y[i] * z[i] * z[i];
    A(i,X2Y2_POLY_IND) = x[i] * x[i] * y[i] * y[i];
    A(i,X2Z2_POLY_IND) = x[i] * x[i] * z[i] * z[i];
    A(i,Y2Z2_POLY_IND) = y[i] * y[i] * z[i] * z[i];
    A(i,X3Y_POLY_IND)  = x[i] * x[i] * x[i] * y[i];
    A(i,X3Z_POLY_IND)  = x[i] * x[i] * x[i] * z[i];
    A(i,Y3X_POLY_IND)  = y[i] * y[i] * y[i] * x[i];
    A(i,Y3Z_POLY_IND)  = y[i] * y[i] * y[i] * z[i];
    A(i,Z3X_POLY_IND)  = z[i] * z[i] * z[i] * x[i];
    A(i,Z3Y_POLY_IND)  = z[i] * z[i] * z[i] * y[i];
    A(i,X4_POLY_IND)   = x[i] * x[i] * x[i] * x[i];
    A(i,Y4_POLY_IND)   = y[i] * y[i] * y[i] * y[i];
    A(i,Z4_POLY_IND)   = z[i] * z[i] * z[i] * z[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(ans(i));
  }
}

DG_FP PolyApprox3D::val_at_2nd(const DG_FP x, const DG_FP y, const DG_FP z) {
  DG_FP res = 0.0;
  res += coeff[C_POLY_IND];
  res += coeff[X_POLY_IND] * x;
  res += coeff[Y_POLY_IND] * y;
  res += coeff[Z_POLY_IND] * z;
  res += coeff[XY_POLY_IND] * x * y;
  res += coeff[XZ_POLY_IND] * x * z;
  res += coeff[YZ_POLY_IND] * y * z;
  res += coeff[X2_POLY_IND] * x * x;
  res += coeff[Y2_POLY_IND] * y * y;
  res += coeff[Z2_POLY_IND] * z * z;
  return res;
}

DG_FP PolyApprox3D::val_at_3rd(const DG_FP x, const DG_FP y, const DG_FP z) {
  DG_FP res = 0.0;
  res += coeff[C_POLY_IND];
  res += coeff[X_POLY_IND] * x;
  res += coeff[Y_POLY_IND] * y;
  res += coeff[Z_POLY_IND] * z;
  res += coeff[XY_POLY_IND] * x * y;
  res += coeff[XZ_POLY_IND] * x * z;
  res += coeff[YZ_POLY_IND] * y * z;
  res += coeff[X2_POLY_IND] * x * x;
  res += coeff[Y2_POLY_IND] * y * y;
  res += coeff[Z2_POLY_IND] * z * z;
  res += coeff[XYZ_POLY_IND] * x * y * z;
  res += coeff[X2Y_POLY_IND] * x * x * y;
  res += coeff[X2Z_POLY_IND] * x * x * z;
  res += coeff[Y2X_POLY_IND] * y * y * x;
  res += coeff[Y2Z_POLY_IND] * y * y * z;
  res += coeff[Z2X_POLY_IND] * z * z * x;
  res += coeff[Z2Y_POLY_IND] * z * z * y;
  res += coeff[X3_POLY_IND] * x * x * x;
  res += coeff[Y3_POLY_IND] * y * y * y;
  res += coeff[Z3_POLY_IND] * z * z * z;
  return res;
}

DG_FP PolyApprox3D::val_at_4th(const DG_FP x, const DG_FP y, const DG_FP z) {
  DG_FP res = 0.0;
  res += coeff[C_POLY_IND];
  res += coeff[X_POLY_IND] * x;
  res += coeff[Y_POLY_IND] * y;
  res += coeff[Z_POLY_IND] * z;
  res += coeff[XY_POLY_IND] * x * y;
  res += coeff[XZ_POLY_IND] * x * z;
  res += coeff[YZ_POLY_IND] * y * z;
  res += coeff[X2_POLY_IND] * x * x;
  res += coeff[Y2_POLY_IND] * y * y;
  res += coeff[Z2_POLY_IND] * z * z;
  res += coeff[XYZ_POLY_IND] * x * y * z;
  res += coeff[X2Y_POLY_IND] * x * x * y;
  res += coeff[X2Z_POLY_IND] * x * x * z;
  res += coeff[Y2X_POLY_IND] * y * y * x;
  res += coeff[Y2Z_POLY_IND] * y * y * z;
  res += coeff[Z2X_POLY_IND] * z * z * x;
  res += coeff[Z2Y_POLY_IND] * z * z * y;
  res += coeff[X3_POLY_IND] * x * x * x;
  res += coeff[Y3_POLY_IND] * y * y * y;
  res += coeff[Z3_POLY_IND] * z * z * z;
  res += coeff[X2YZ_POLY_IND] * x * x * y * z;
  res += coeff[XY2Z_POLY_IND] * x * y * y * z;
  res += coeff[XYZ2_POLY_IND] * x * y * z * z;
  res += coeff[X2Y2_POLY_IND] * x * x * y * y;
  res += coeff[X2Z2_POLY_IND] * x * x * z * z;
  res += coeff[Y2Z2_POLY_IND] * y * y * z * z;
  res += coeff[X3Y_POLY_IND] * x * x * x * y;
  res += coeff[X3Z_POLY_IND] * x * x * x * z;
  res += coeff[Y3X_POLY_IND] * y * y * y * x;
  res += coeff[Y3Z_POLY_IND] * y * y * y * z;
  res += coeff[Z3X_POLY_IND] * z * z * z * x;
  res += coeff[Z3Y_POLY_IND] * z * z * z * y;
  res += coeff[X4_POLY_IND] * x * x * x * x;
  res += coeff[Y4_POLY_IND] * y * y * y * y;
  res += coeff[Z4_POLY_IND] * z * z * z * z;
  return res;
}

DG_FP PolyApprox3D::val_at(const DG_FP x, const DG_FP y, const DG_FP z) {
  DG_FP res = 0.0;
  if(N == 2) {
    res = val_at_2nd(x - offset_x, y - offset_y, z - offset_z);
  } else if(N == 3) {
    res = val_at_3rd(x - offset_x, y - offset_y, z - offset_z);
  } else if(N == 4) {
    res = val_at_4th(x - offset_x, y - offset_y, z - offset_z);
  }
  return res;
}

void PolyApprox3D::grad_at_2nd(const DG_FP x, const DG_FP y, const DG_FP z,
                             DG_FP &dx, DG_FP &dy, DG_FP &dz) {
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;

  dx += coeff[X_POLY_IND];
  dx += coeff[XY_POLY_IND] * y;
  dx += coeff[XZ_POLY_IND] * z;
  dx += 2.0 * coeff[X2_POLY_IND] * x;

  dy += coeff[Y_POLY_IND];
  dy += coeff[XY_POLY_IND] * x;
  dy += coeff[YZ_POLY_IND] * z;
  dy += 2.0 * coeff[Y2_POLY_IND] * y;

  dz += coeff[Z_POLY_IND];
  dz += coeff[XZ_POLY_IND] * x;
  dz += coeff[YZ_POLY_IND] * y;
  dz += 2.0 * coeff[Z2_POLY_IND] * z;
}

void PolyApprox3D::grad_at_3rd(const DG_FP x, const DG_FP y, const DG_FP z,
                             DG_FP &dx, DG_FP &dy, DG_FP &dz) {
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;

  dx += coeff[X_POLY_IND];
  dx += coeff[XY_POLY_IND] * y;
  dx += coeff[XZ_POLY_IND] * z;
  dx += 2.0 * coeff[X2_POLY_IND] * x;
  dx += coeff[XYZ_POLY_IND] * y * z;
  dx += 2.0 * coeff[X2Y_POLY_IND] * x * y;
  dx += 2.0 * coeff[X2Z_POLY_IND] * x * z;
  dx += coeff[Y2X_POLY_IND] * y * y;
  dx += coeff[Z2X_POLY_IND] * z * z;
  dx += 3.0 * coeff[X3_POLY_IND] * x * x;

  dy += coeff[Y_POLY_IND];
  dy += coeff[XY_POLY_IND] * x;
  dy += coeff[YZ_POLY_IND] * z;
  dy += 2.0 * coeff[Y2_POLY_IND] * y;
  dy += coeff[XYZ_POLY_IND] * x * z;
  dy += coeff[X2Y_POLY_IND] * x * x;
  dy += 2.0 * coeff[Y2X_POLY_IND] * y * x;
  dy += 2.0 * coeff[Y2Z_POLY_IND] * y * z;
  dy += coeff[Z2Y_POLY_IND] * z * z;
  dy += 3.0 * coeff[Y3_POLY_IND] * y * y;

  dz += coeff[Z_POLY_IND];
  dz += coeff[XZ_POLY_IND] * x;
  dz += coeff[YZ_POLY_IND] * y;
  dz += 2.0 * coeff[Z2_POLY_IND] * z;
  dz += coeff[XYZ_POLY_IND] * x * y;
  dz += coeff[X2Z_POLY_IND] * x * x;
  dz += coeff[Y2Z_POLY_IND] * y * y;
  dz += 2.0 * coeff[Z2X_POLY_IND] * z * x;
  dz += 2.0 * coeff[Z2Y_POLY_IND] * z * y;
  dz += 3.0 * coeff[Z3_POLY_IND] * z * z;
}

void PolyApprox3D::grad_at_4th(const DG_FP x, const DG_FP y, const DG_FP z,
                             DG_FP &dx, DG_FP &dy, DG_FP &dz) {
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;

  dx += coeff[X_POLY_IND];
  dx += coeff[XY_POLY_IND] * y;
  dx += coeff[XZ_POLY_IND] * z;
  dx += 2.0 * coeff[X2_POLY_IND] * x;
  dx += coeff[XYZ_POLY_IND] * y * z;
  dx += 2.0 * coeff[X2Y_POLY_IND] * x * y;
  dx += 2.0 * coeff[X2Z_POLY_IND] * x * z;
  dx += coeff[Y2X_POLY_IND] * y * y;
  dx += coeff[Z2X_POLY_IND] * z * z;
  dx += 3.0 * coeff[X3_POLY_IND] * x * x;
  dx += 2.0 * coeff[X2YZ_POLY_IND] * x * y * z;
  dx += coeff[XY2Z_POLY_IND] * y * y * z;
  dx += coeff[XYZ2_POLY_IND] * y * z * z;
  dx += 2.0 * coeff[X2Y2_POLY_IND] * x * y * y;
  dx += 2.0 * coeff[X2Z2_POLY_IND] * x * z * z;
  dx += 3.0 * coeff[X3Y_POLY_IND] * x * x * y;
  dx += 3.0 * coeff[X3Z_POLY_IND] * x * x * z;
  dx += coeff[Y3X_POLY_IND] * y * y * y;
  dx += coeff[Z3X_POLY_IND] * z * z * z;
  dx += 4.0 * coeff[X4_POLY_IND] * x * x * x;

  dy += coeff[Y_POLY_IND];
  dy += coeff[XY_POLY_IND] * x;
  dy += coeff[YZ_POLY_IND] * z;
  dy += 2.0 * coeff[Y2_POLY_IND] * y;
  dy += coeff[XYZ_POLY_IND] * x * z;
  dy += coeff[X2Y_POLY_IND] * x * x;
  dy += 2.0 * coeff[Y2X_POLY_IND] * y * x;
  dy += 2.0 * coeff[Y2Z_POLY_IND] * y * z;
  dy += coeff[Z2Y_POLY_IND] * z * z;
  dy += 3.0 * coeff[Y3_POLY_IND] * y * y;
  dy += coeff[X2YZ_POLY_IND] * x * x * z;
  dy += 2.0 * coeff[XY2Z_POLY_IND] * x * y * z;
  dy += coeff[XYZ2_POLY_IND] * x * z * z;
  dy += 2.0 * coeff[X2Y2_POLY_IND] * x * x * y;
  dy += 2.0 * coeff[Y2Z2_POLY_IND] * y * z * z;
  dy += coeff[X3Y_POLY_IND] * x * x * x;
  dy += 3.0 * coeff[Y3X_POLY_IND] * y * y * x;
  dy += 3.0 * coeff[Y3Z_POLY_IND] * y * y * z;
  dy += coeff[Z3Y_POLY_IND] * z * z * z;
  dy += 4.0 * coeff[Y4_POLY_IND] * y * y * y;

  dz += coeff[Z_POLY_IND];
  dz += coeff[XZ_POLY_IND] * x;
  dz += coeff[YZ_POLY_IND] * y;
  dz += 2.0 * coeff[Z2_POLY_IND] * z;
  dz += coeff[XYZ_POLY_IND] * x * y;
  dz += coeff[X2Z_POLY_IND] * x * x;
  dz += coeff[Y2Z_POLY_IND] * y * y;
  dz += 2.0 * coeff[Z2X_POLY_IND] * z * x;
  dz += 2.0 * coeff[Z2Y_POLY_IND] * z * y;
  dz += 3.0 * coeff[Z3_POLY_IND] * z * z;
  dz += coeff[X2YZ_POLY_IND] * x * x * y;
  dz += coeff[XY2Z_POLY_IND] * x * y * y;
  dz += 2.0 * coeff[XYZ2_POLY_IND] * x * y * z;
  dz += 2.0 * coeff[X2Z2_POLY_IND] * x * x * z;
  dz += 2.0 * coeff[Y2Z2_POLY_IND] * y * y * z;
  dz += coeff[X3Z_POLY_IND] * x * x * x;
  dz += coeff[Y3Z_POLY_IND] * y * y * y;
  dz += 3.0 * coeff[Z3X_POLY_IND] * z * z * x;
  dz += 3.0 * coeff[Z3Y_POLY_IND] * z * z * y;
  dz += 4.0 * coeff[Z4_POLY_IND] * z * z * z;
}

void PolyApprox3D::grad_at(const DG_FP x, const DG_FP y, const DG_FP z,
                         DG_FP &dx, DG_FP &dy, DG_FP &dz) {
  if(N == 2) {
    grad_at_2nd(x - offset_x, y - offset_y, z - offset_z, dx, dy, dz);
  } else if(N == 3) {
    grad_at_3rd(x - offset_x, y - offset_y, z - offset_z, dx, dy, dz);
  } else if(N == 4) {
    grad_at_4th(x - offset_x, y - offset_y, z - offset_z, dx, dy, dz);
  }
}

void PolyApprox3D::hessian_at_2nd(const DG_FP x, const DG_FP y, const DG_FP z,
                                DG_FP &dx2, DG_FP &dy2, DG_FP &dz2,
                                DG_FP &dxy, DG_FP &dxz, DG_FP &dyz) {
  dx2 = 0.0;
  dy2 = 0.0;
  dz2 = 0.0;
  dxy = 0.0;
  dxz = 0.0;
  dyz = 0.0;

  dx2 += 2.0 * coeff[X2_POLY_IND];
  dy2 += 2.0 * coeff[Y2_POLY_IND];
  dz2 += 2.0 * coeff[Z2_POLY_IND];
  dxy += coeff[XY_POLY_IND];
  dxz += coeff[XZ_POLY_IND];
  dyz += coeff[YZ_POLY_IND];
}

void PolyApprox3D::hessian_at_3rd(const DG_FP x, const DG_FP y, const DG_FP z,
                                DG_FP &dx2, DG_FP &dy2, DG_FP &dz2,
                                DG_FP &dxy, DG_FP &dxz, DG_FP &dyz) {
  dx2 = 0.0;
  dy2 = 0.0;
  dz2 = 0.0;
  dxy = 0.0;
  dxz = 0.0;
  dyz = 0.0;

  dx2 += 2.0 * coeff[X2_POLY_IND];
  dx2 += 2.0 * coeff[X2Y_POLY_IND] * y;
  dx2 += 2.0 * coeff[X2Z_POLY_IND] * z;
  dx2 += 6.0 * coeff[X3_POLY_IND] * x;

  dy2 += 2.0 * coeff[Y2_POLY_IND];
  dy2 += 2.0 * coeff[Y2X_POLY_IND] * x;
  dy2 += 2.0 * coeff[Y2Z_POLY_IND] * z;
  dy2 += 6.0 * coeff[Y3_POLY_IND] * y;

  dz2 += 2.0 * coeff[Z2_POLY_IND];
  dz2 += 2.0 * coeff[Z2X_POLY_IND] * x;
  dz2 += 2.0 * coeff[Z2Y_POLY_IND] * y;
  dz2 += 6.0 * coeff[Z3_POLY_IND] * z;

  dxy += coeff[XY_POLY_IND];
  dxy += coeff[XYZ_POLY_IND] * z;
  dxy += 2.0 * coeff[X2Y_POLY_IND] * x;
  dxy += 2.0 * coeff[Y2X_POLY_IND] * y;

  dxz += coeff[XZ_POLY_IND];
  dxz += coeff[XYZ_POLY_IND] * y;
  dxz += 2.0 * coeff[X2Z_POLY_IND] * x;
  dxz += 2.0 * coeff[Z2X_POLY_IND] * z;

  dyz += coeff[YZ_POLY_IND];
  dyz += coeff[XYZ_POLY_IND] * x;
  dyz += 2.0 * coeff[Y2Z_POLY_IND] * y;
  dyz += 2.0 * coeff[Z2Y_POLY_IND] * z;
}

void PolyApprox3D::hessian_at_4th(const DG_FP x, const DG_FP y, const DG_FP z,
                                DG_FP &dx2, DG_FP &dy2, DG_FP &dz2,
                                DG_FP &dxy, DG_FP &dxz, DG_FP &dyz) {
  dx2 = 0.0;
  dy2 = 0.0;
  dz2 = 0.0;
  dxy = 0.0;
  dxz = 0.0;
  dyz = 0.0;

  dx2 += 2.0 * coeff[X2_POLY_IND];
  dx2 += 2.0 * coeff[X2Y_POLY_IND] * y;
  dx2 += 2.0 * coeff[X2Z_POLY_IND] * z;
  dx2 += 6.0 * coeff[X3_POLY_IND] * x;
  dx2 += 2.0 * coeff[X2YZ_POLY_IND] * y * z;
  dx2 += 2.0 * coeff[X2Y2_POLY_IND] * y * y;
  dx2 += 2.0 * coeff[X2Z2_POLY_IND] * z * z;
  dx2 += 6.0 * coeff[X3Y_POLY_IND] * x * y;
  dx2 += 6.0 * coeff[X3Z_POLY_IND] * x * z;
  dx2 += 12.0 * coeff[X4_POLY_IND] * x * x;

  dy2 += 2.0 * coeff[Y2_POLY_IND];
  dy2 += 2.0 * coeff[Y2X_POLY_IND] * x;
  dy2 += 2.0 * coeff[Y2Z_POLY_IND] * z;
  dy2 += 6.0 * coeff[Y3_POLY_IND] * y;
  dy2 += 2.0 * coeff[XY2Z_POLY_IND] * x * z;
  dy2 += 2.0 * coeff[X2Y2_POLY_IND] * x * x;
  dy2 += 2.0 * coeff[Y2Z2_POLY_IND] * z * z;
  dy2 += 6.0 * coeff[Y3X_POLY_IND] * y * x;
  dy2 += 6.0 * coeff[Y3Z_POLY_IND] * y * z;
  dy2 += 12.0 * coeff[Y4_POLY_IND] * y * y;

  dz2 += 2.0 * coeff[Z2_POLY_IND];
  dz2 += 2.0 * coeff[Z2X_POLY_IND] * x;
  dz2 += 2.0 * coeff[Z2Y_POLY_IND] * y;
  dz2 += 6.0 * coeff[Z3_POLY_IND] * z;
  dz2 += 2.0 * coeff[XYZ2_POLY_IND] * x * y;
  dz2 += 2.0 * coeff[X2Z2_POLY_IND] * x * x;
  dz2 += 2.0 * coeff[Y2Z2_POLY_IND] * y * y;
  dz2 += 6.0 * coeff[Z3X_POLY_IND] * z * x;
  dz2 += 6.0 * coeff[Z3Y_POLY_IND] * z * y;
  dz2 += 12.0 * coeff[Z4_POLY_IND] * z * z;

  dxy += coeff[XY_POLY_IND];
  dxy += coeff[XYZ_POLY_IND] * z;
  dxy += 2.0 * coeff[X2Y_POLY_IND] * x;
  dxy += 2.0 * coeff[Y2X_POLY_IND] * y;
  dxy += 2.0 * coeff[X2YZ_POLY_IND] * x * z;
  dxy += 2.0 * coeff[XY2Z_POLY_IND] * y * z;
  dxy += coeff[XYZ2_POLY_IND] * z * z;
  dxy += 4.0 * coeff[X2Y2_POLY_IND] * x * y;
  dxy += 3.0 * coeff[X3Y_POLY_IND] * x * x;
  dxy += 3.0 * coeff[Y3X_POLY_IND] * y * y;

  dxz += coeff[XZ_POLY_IND];
  dxz += coeff[XYZ_POLY_IND] * y;
  dxz += 2.0 * coeff[X2Z_POLY_IND] * x;
  dxz += 2.0 * coeff[Z2X_POLY_IND] * z;
  dxz += 2.0 * coeff[X2YZ_POLY_IND] * x * y;
  dxz += coeff[XY2Z_POLY_IND] * y * y;
  dxz += 2.0 * coeff[XYZ2_POLY_IND] * y * z;
  dxz += 4.0 * coeff[X2Z2_POLY_IND] * x * z;
  dxz += 3.0 * coeff[X3Z_POLY_IND] * x * x;
  dxz += 3.0 * coeff[Z3X_POLY_IND] * z * z;

  dyz += coeff[YZ_POLY_IND];
  dyz += coeff[XYZ_POLY_IND] * x;
  dyz += 2.0 * coeff[Y2Z_POLY_IND] * y;
  dyz += 2.0 * coeff[Z2Y_POLY_IND] * z;
  dyz += coeff[X2YZ_POLY_IND] * x * x;
  dyz += 2.0 * coeff[XY2Z_POLY_IND] * x * y;
  dyz += 2.0 * coeff[XYZ2_POLY_IND] * x * z;
  dyz += 4.0 * coeff[Y2Z2_POLY_IND] * y * z;
  dyz += 3.0 * coeff[Y3Z_POLY_IND] * y * y;
  dyz += 3.0 * coeff[Z3Y_POLY_IND] * z * z;
}

void PolyApprox3D::hessian_at(const DG_FP x, const DG_FP y, const DG_FP z,
                            DG_FP &dx2, DG_FP &dy2, DG_FP &dz2,
                            DG_FP &dxy, DG_FP &dxz, DG_FP &dyz) {
  if(N == 2) {
    hessian_at_2nd(x - offset_x, y - offset_y, z - offset_z, dx2, dy2, dz2, dxy, dxz, dyz);
  } else if(N == 3) {
    hessian_at_3rd(x - offset_x, y - offset_y, z - offset_z, dx2, dy2, dz2, dxy, dxz, dyz);
  } else if(N == 4) {
    hessian_at_4th(x - offset_x, y - offset_y, z - offset_z, dx2, dy2, dz2, dxy, dxz, dyz);
  }
}

int PolyApprox3D::num_coeff() {
  if(N == 2) {
    return 10;
  } else if(N == 3) {
    return 20;
  } else if(N == 4) {
    return 35;
  } else {
    return -1;
  }
}

int PolyApprox3D::num_pts() {
  if(N == 2) {
    return 10;
  } else if(N == 3) {
    return 20;
  } else if(N == 4) {
    return 35;
  } else {
    return 0;
  }
}

int PolyApprox3D::num_elem_stencil() {
  if(N == 2) {
    return 12;
  } else if(N == 3) {
    return 17;
  } else if(N == 4) {
    return 8;
  } else {
    return 0;
  }
}

DG_FP PolyApprox3D::get_coeff(int ind) {
  return coeff[ind];
}

void PolyApprox3D::get_offsets(DG_FP &x, DG_FP &y, DG_FP &z) {
  x = offset_x;
  y = offset_y;
  z = offset_z;
}

struct stencil_query {
  int ind;
  set<int> central_inds;
};

map<int,set<int>> PolyApprox3D::get_stencils(const set<int> &central_inds, op_map edge_map) {
  timer->startTimer("PolyApprox3D - get_stencils");
  map<int,set<int>> stencils;
  map<int,stencil_query> queryInds;
  const int num_elements = num_elem_stencil();

  for(const auto &ind : central_inds) {
    set<int> st;
    st.insert(ind);
    stencils.insert({ind, st});
    stencil_query sq;
    sq.ind = ind;
    sq.central_inds.insert(ind);
    queryInds.insert({ind, sq});
  }

  if(num_elements > 0) {
    const int numEdges = edge_map->from->size;
    while(queryInds.size() > 0) {
      map<int,stencil_query> newQueryInds;

      // Iterate over each edge pair
      for(int i = 0; i < numEdges * 2; i++) {
        // Find if this cell ind is in the query inds
        auto it = queryInds.find(edge_map->map[i]);
        if(it != queryInds.end()) {
          if(i % 2 == 0) {
            // For each central ind associated with this query ind
            for(const auto &ind : it->second.central_inds) {
              auto stencil_it = stencils.find(ind);
              // Check if the other cell in this edge is already in the stencil for this central ind
              if(stencil_it->second.find(edge_map->map[i + 1]) == stencil_it->second.end()
                 && stencil_it->second.size() < num_elements) {
                stencil_it->second.insert(edge_map->map[i + 1]);
                // If stencil is not full then add to next rounds query inds
                if(stencil_it->second.size() < num_elements) {
                  stencil_query sq;
                  sq.ind = edge_map->map[i + 1];
                  auto res = newQueryInds.insert({edge_map->map[i + 1], sq});
                  res.first->second.central_inds.insert(ind);
                }
              }
            }
          } else {
            // For each central ind associated with this query ind
            for(const auto &ind : it->second.central_inds) {
              auto stencil_it = stencils.find(ind);
              // Check if the other cell in this edge is already in the stencil for this central ind
              if(stencil_it->second.find(edge_map->map[i - 1]) == stencil_it->second.end()
                 && stencil_it->second.size() < num_elements) {
                stencil_it->second.insert(edge_map->map[i - 1]);
                // If stencil is not full then add to next rounds query inds
                if(stencil_it->second.size() < num_elements) {
                  stencil_query sq;
                  sq.ind = edge_map->map[i - 1];
                  auto res = newQueryInds.insert({edge_map->map[i - 1], sq});
                  res.first->second.central_inds.insert(ind);
                }
              }
            }
          }
        }
      }

      queryInds = newQueryInds;
    }
  }
  timer->endTimer("PolyApprox3D - get_stencils");
  return stencils;
}
