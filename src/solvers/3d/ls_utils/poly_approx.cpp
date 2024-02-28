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
                           const DG_FP *z_ptr, const DG_FP *s_ptr, 
                           const DG_FP h) {
  get_offset(cell_ind, x_ptr, y_ptr, z_ptr);

  vector<DG_FP> x_vec, y_vec, z_vec, s_vec;
  stencil_data(cell_ind, stencil, x_ptr, y_ptr, z_ptr, s_ptr, x_vec, y_vec, z_vec, s_vec);

  fit_poly(x_vec, y_vec, z_vec, s_vec, h);
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
  offset_x = 0.0;
  offset_y = 0.0;
  offset_z = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    offset_x += x_ptr[ind * DG_NP + i];
    offset_y += y_ptr[ind * DG_NP + i];
    offset_z += z_ptr[ind * DG_NP + i];
  }
  offset_x /= (DG_FP)DG_NP;
  offset_y /= (DG_FP)DG_NP;
  offset_z /= (DG_FP)DG_NP;
}

void PolyApprox::fit_poly(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &z, const vector<DG_FP> &s, const DG_FP h) {
  // Set A vandermonde matrix and b
  arma::mat A = get_vandermonde(x, y, z);
  arma::vec b(s);

  // Fit poly by using QR
  arma::mat Q, R;
  arma::qr(Q, R, A);
  arma::vec ans = arma::solve(R, Q.t() * b);
  
  coeff = arma::conv_to<vector<DG_FP>>::from(ans);

  arma::vec res = A * ans;
  bool node_with_wrong_sign = false;
  set<int> problem_inds;
  arma::vec w(x.size());
  const DG_FP sigma = h * 0.01;
  for(int i = 0; i < x.size(); i++) {
    w(i) = exp(-s[i] * s[i] / sigma);
    if(res(i) > 0.0 != b(i) > 0.0) {
      node_with_wrong_sign = true;
      problem_inds.insert(i);
    }
  }

  int redo_counter = 0;
  const int max_redo = 50;
  while(node_with_wrong_sign && redo_counter < max_redo) {
    for(int i = 0; i < x.size(); i++) {
      if(problem_inds.count(i) > 0) 
        w(i) = w(i) * 2.0;
    }
    arma::mat W = arma::diagmat(arma::sqrt(w));
    arma::mat Q, R;
    arma::qr(Q, R, W*A);
    arma::vec ans = arma::solve(R, Q.t() * W * b);
  
    coeff = arma::conv_to<vector<DG_FP>>::from(ans);

    arma::vec res = A * ans;
    node_with_wrong_sign = false;
    problem_inds.clear();
    for(int i = 0; i < x.size(); i++) {
      if(res(i) > 0.0 != b(i) > 0.0) {
        node_with_wrong_sign = true;
        problem_inds.insert(i);
      }
    }

    redo_counter++;
  }

  // if(redo_counter == max_redo) printf("Max redo\n");
  // if(node_with_wrong_sign) printf("Node with wrong sign\n");
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
  // Setup random number generator for later
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
        // if(sten == cell_ind) {
        //   res.first->second.val = s_ptr[ind];
        // }
        // Point already exists
        res.first->second.val += s_ptr[ind];
        res.first->second.count++;
      }
    }
  }

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

arma::mat PolyApprox3D::get_2nd_order_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y, 
                                                  const vector<DG_FP> &z) {
  arma::mat A(x.size(), num_coeff());
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
  }
  return A;
}

arma::mat PolyApprox3D::get_3rd_order_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y,
                                                  const vector<DG_FP> &z) {
  arma::mat A(x.size(), num_coeff());
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
  }
  return A;
}

arma::mat PolyApprox3D::get_4th_order_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y,
                                                  const vector<DG_FP> &z) {
  arma::mat A(x.size(), num_coeff());
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
  }
  return A;
}

arma::mat PolyApprox::get_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &z) {
  arma::mat res;
  if(N == 2) {
    res = get_2nd_order_vandermonde(x, y, z);
  } else if(N == 3) {
    res = get_3rd_order_vandermonde(x, y, z);
  } else if(N == 4) {
    res = get_4th_order_vandermonde(x, y, z);
  }
  return res;
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
    res = val_at_2nd(x, y, z);
  } else if(N == 3) {
    res = val_at_3rd(x, y, z);
  } else if(N == 4) {
    res = val_at_4th(x, y, z);
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
    grad_at_2nd(x, y, z, dx, dy, dz);
  } else if(N == 3) {
    grad_at_3rd(x, y, z, dx, dy, dz);
  } else if(N == 4) {
    grad_at_4th(x, y, z, dx, dy, dz);
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
    hessian_at_2nd(x, y, z, dx2, dy2, dz2, dxy, dxz, dyz);
  } else if(N == 3) {
    hessian_at_3rd(x, y, z, dx2, dy2, dz2, dxy, dxz, dyz);
  } else if(N == 4) {
    hessian_at_4th(x, y, z, dx2, dy2, dz2, dxy, dxz, dyz);
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
    return 30;
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
