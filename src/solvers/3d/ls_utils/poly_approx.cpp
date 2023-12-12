#include "ls_utils/3d/poly_approx.h"

#include "timing.h"

extern Timing *timer;

using namespace std;

void num_pts_pos_neg(const vector<DG_FP> &s, int &pos, int &neg);

PolyApprox3D::PolyApprox3D(const int cell_ind, set<int> stencil,
                       const DG_FP *x_ptr, const DG_FP *y_ptr,
                       const DG_FP *z_ptr, const DG_FP *s_ptr) {
  get_offset(cell_ind, x_ptr, y_ptr, z_ptr);

  vector<DG_FP> x_vec, y_vec, z_vec, s_vec;
  stencil_data(cell_ind, stencil, x_ptr, y_ptr, z_ptr, s_ptr, x_vec, y_vec, z_vec, s_vec);

  // Make sure equal number of points on each side of the line
/*
  int pts_needed = num_coeff();
  int pts_aim = num_pts();
  int pos_pts, neg_pts;
  num_pts_pos_neg(s_vec, pos_pts, neg_pts);
  while(x_vec.size() > pts_needed) {
    // Find point furthest from the interface to discard
    int ind_discard;
    if(pos_pts > neg_pts) {
      DG_FP max = s_vec[0];
      ind_discard = 0;
      for(int i = 1; i < x_vec.size(); i++) {
        if(s_vec[i] > max) {
          max = s_vec[i];
          ind_discard = i;
        }
      }
    } else {
      DG_FP min = s_vec[0];
      ind_discard = 0;
      for(int i = 1; i < x_vec.size(); i++) {
        if(s_vec[i] < min) {
          min = s_vec[i];
          ind_discard = i;
        }
      }
    }
    // Discard ind
    x_vec.erase(x_vec.begin() + ind_discard);
    y_vec.erase(y_vec.begin() + ind_discard);
    z_vec.erase(z_vec.begin() + ind_discard);
    s_vec.erase(s_vec.begin() + ind_discard);

    num_pts_pos_neg(s_vec, pos_pts, neg_pts);
  }
*/

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
  // offset_x = x_ptr[ind * DG_NP];
  // offset_y = y_ptr[ind * DG_NP];
  // offset_z = z_ptr[ind * DG_NP];
  offset_x = 0.0;
  offset_y = 0.0;
  offset_z = 0.0;
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

void PolyApprox3D::stencil_data(const int cell_ind, const set<int> &stencil,
                              const DG_FP *x_ptr, const DG_FP *y_ptr,
                              const DG_FP *z_ptr, const DG_FP *s_ptr,
                              vector<DG_FP> &x, vector<DG_FP> &y,
                              vector<DG_FP> &z, vector<DG_FP> &s) {
  map<Coord, Point, cmpCoords> pointMap;

  for(const auto &sten : stencil) {
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
        if(sten == cell_ind) {
          res.first->second.val = s_ptr[ind];
        }
        // Point already exists
        // res.first->second.val += s_ptr[ind];
        // res.first->second.count++;
      }
    }
  }

  for(auto const &p : pointMap) {
    x.push_back(p.second.coord.x);
    y.push_back(p.second.coord.y);
    z.push_back(p.second.coord.z);
    s.push_back(p.second.val);
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

void num_pts_pos_neg(const vector<DG_FP> &s, int &pos, int &neg) {
  pos = 0;
  neg = 0;
  for(int i = 0; i < s.size(); i++) {
    if(s[i] > 0.0) pos++;
    if(s[i] < 0.0) neg++;
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
    return 0;
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
