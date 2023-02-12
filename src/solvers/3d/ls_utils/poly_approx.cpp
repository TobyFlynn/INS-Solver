#include "ls_utils/3d/poly_approx.h"

#include "timing.h"

extern Timing *timer;

using namespace std;

void num_pts_pos_neg(const vector<double> &s, int &pos, int &neg);

PolyApprox3D::PolyApprox3D(const int cell_ind, set<int> stencil,
                       const double *x_ptr, const double *y_ptr,
                       const double *z_ptr, const double *s_ptr) {
  get_offset(cell_ind, x_ptr, y_ptr, z_ptr);

  vector<double> x_vec, y_vec, z_vec, s_vec;
  stencil_data(cell_ind, stencil, x_ptr, y_ptr, z_ptr, s_ptr, x_vec, y_vec, z_vec, s_vec);

  // Make sure equal number of points on each side of the line
  // int pts_needed = num_coeff();
  // int pts_aim = num_pts();
  // int pos_pts, neg_pts;
  // num_pts_pos_neg(s_vec, pos_pts, neg_pts);
  // while((x_vec.size() > pts_needed && pos_pts != neg_pts) || x_vec.size() > pts_aim) {
  //   // Find point furthest from the interface to discard
  //   int ind_discard;
  //   if(pos_pts > neg_pts) {
  //     double max = s_vec[0];
  //     ind_discard = 0;
  //     for(int i = 1; i < x_vec.size(); i++) {
  //       if(s_vec[i] > max) {
  //         max = s_vec[i];
  //         ind_discard = i;
  //       }
  //     }
  //   } else {
  //     double min = s_vec[0];
  //     ind_discard = 0;
  //     for(int i = 1; i < x_vec.size(); i++) {
  //       if(s_vec[i] < min) {
  //         min = s_vec[i];
  //         ind_discard = i;
  //       }
  //     }
  //   }
  //   // Discard ind
  //   x_vec.erase(x_vec.begin() + ind_discard);
  //   y_vec.erase(y_vec.begin() + ind_discard);
  //   z_vec.erase(z_vec.begin() + ind_discard);
  //   s_vec.erase(s_vec.begin() + ind_discard);
  //
  //   num_pts_pos_neg(s_vec, pos_pts, neg_pts);
  // }

  if(N == 2) {
    set_2nd_order_coeff(x_vec, y_vec, z_vec, s_vec);
  } else if(N == 3) {
    set_3rd_order_coeff(x_vec, y_vec, z_vec, s_vec);
  } else if(N == 4) {
    set_4th_order_coeff(x_vec, y_vec, z_vec, s_vec);
  }
}

PolyApprox3D::PolyApprox3D(std::vector<double> &c, double off_x, double off_y,
                       double off_z) {
  offset_x = off_x;
  offset_y = off_y;
  offset_z = off_z;
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(c[i]);
  }
}

void PolyApprox3D::get_offset(const int ind, const double *x_ptr, const double *y_ptr,
                            const double *z_ptr) {
  // offset_x = x_ptr[ind * DG_NP];
  // offset_y = y_ptr[ind * DG_NP];
  // offset_z = z_ptr[ind * DG_NP];
  offset_x = 0.0;
  offset_y = 0.0;
  offset_z = 0.0;
}

struct Coord {
  double x;
  double y;
  double z;
};

struct Point {
  Coord coord;
  double val;
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
                              const double *x_ptr, const double *y_ptr,
                              const double *z_ptr, const double *s_ptr,
                              vector<double> &x, vector<double> &y,
                              vector<double> &z, vector<double> &s) {
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
    s.push_back(p.second.val / (double)p.second.count);
  }
}

void PolyApprox3D::set_2nd_order_coeff(const vector<double> &x, const vector<double> &y,
                                     const vector<double> &z, const vector<double> &s) {
  arma::mat A(x.size(), num_coeff());
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,0) = 1.0;
    A(i,1) = x[i];
    A(i,2) = y[i];
    A(i,3) = z[i];
    A(i,4) = x[i] * y[i];
    A(i,5) = x[i] * z[i];
    A(i,6) = y[i] * z[i];
    A(i,7) = x[i] * x[i];
    A(i,8) = y[i] * y[i];
    A(i,9) = z[i] * z[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(ans(i));
  }
}

void PolyApprox3D::set_3rd_order_coeff(const vector<double> &x, const vector<double> &y,
                                     const vector<double> &z, const vector<double> &s) {
  arma::mat A(x.size(), num_coeff());
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,0) = 1.0;
    A(i,1) = x[i];
    A(i,2) = y[i];
    A(i,3) = z[i];
    A(i,4) = x[i] * y[i];
    A(i,5) = x[i] * z[i];
    A(i,6) = y[i] * z[i];
    A(i,7) = x[i] * x[i];
    A(i,8) = y[i] * y[i];
    A(i,9) = z[i] * z[i];
    A(i,10) = x[i] * x[i] * y[i];
    A(i,11) = x[i] * x[i] * z[i];
    A(i,12) = y[i] * y[i] * x[i];
    A(i,13) = y[i] * y[i] * z[i];
    A(i,14) = z[i] * z[i] * x[i];
    A(i,15) = z[i] * z[i] * y[i];
    A(i,16) = x[i] * x[i] * x[i];
    A(i,17) = y[i] * y[i] * y[i];
    A(i,18) = z[i] * z[i] * z[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(ans(i));
  }
}

void PolyApprox3D::set_4th_order_coeff(const vector<double> &x, const vector<double> &y,
                                     const vector<double> &z, const vector<double> &s) {
  arma::mat A(x.size(), num_coeff());
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,0) = 1.0;
    A(i,1) = x[i];
    A(i,2) = y[i];
    A(i,3) = z[i];
    A(i,4) = x[i] * y[i];
    A(i,5) = x[i] * z[i];
    A(i,6) = y[i] * z[i];
    A(i,7) = x[i] * x[i];
    A(i,8) = y[i] * y[i];
    A(i,9) = z[i] * z[i];
    A(i,10) = x[i] * x[i] * y[i];
    A(i,11) = x[i] * x[i] * z[i];
    A(i,12) = y[i] * y[i] * x[i];
    A(i,13) = y[i] * y[i] * z[i];
    A(i,14) = z[i] * z[i] * x[i];
    A(i,15) = z[i] * z[i] * y[i];
    A(i,16) = x[i] * x[i] * x[i];
    A(i,17) = y[i] * y[i] * y[i];
    A(i,18) = z[i] * z[i] * z[i];
    A(i,19) = x[i] * x[i] * y[i] * y[i];
    A(i,20) = x[i] * x[i] * z[i] * z[i];
    A(i,21) = y[i] * y[i] * z[i] * z[i];
    A(i,22) = x[i] * x[i] * x[i] * y[i];
    A(i,23) = x[i] * x[i] * x[i] * z[i];
    A(i,24) = y[i] * y[i] * y[i] * x[i];
    A(i,25) = y[i] * y[i] * y[i] * z[i];
    A(i,26) = z[i] * z[i] * z[i] * x[i];
    A(i,27) = z[i] * z[i] * z[i] * y[i];
    A(i,28) = x[i] * x[i] * x[i] * x[i];
    A(i,29) = y[i] * y[i] * y[i] * y[i];
    A(i,30) = z[i] * z[i] * z[i] * z[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(ans(i));
  }
}

double PolyApprox3D::val_at_2nd(const double x, const double y, const double z) {
  double res = 0.0;
  res += coeff[0];
  res += coeff[1] * x;
  res += coeff[2] * y;
  res += coeff[3] * z;
  res += coeff[4] * x * y;
  res += coeff[5] * x * z;
  res += coeff[6] * y * z;
  res += coeff[7] * x * x;
  res += coeff[8] * y * y;
  res += coeff[9] * z * z;
  return res;
}

double PolyApprox3D::val_at_3rd(const double x, const double y, const double z) {
  double res = 0.0;
  res += coeff[0];
  res += coeff[1] * x;
  res += coeff[2] * y;
  res += coeff[3] * z;
  res += coeff[4] * x * y;
  res += coeff[5] * x * z;
  res += coeff[6] * y * z;
  res += coeff[7] * x * x;
  res += coeff[8] * y * y;
  res += coeff[9] * z * z;
  res += coeff[10] * x * x * y;
  res += coeff[11] * x * x * z;
  res += coeff[12] * y * y * x;
  res += coeff[13] * y * y * z;
  res += coeff[14] * z * z * x;
  res += coeff[15] * z * z * y;
  res += coeff[16] * x * x * x;
  res += coeff[17] * y * y * y;
  res += coeff[18] * z * z * z;
  return res;
}

double PolyApprox3D::val_at_4th(const double x, const double y, const double z) {
  double res = 0.0;
  res += coeff[0];
  res += coeff[1] * x;
  res += coeff[2] * y;
  res += coeff[3] * z;
  res += coeff[4] * x * y;
  res += coeff[5] * x * z;
  res += coeff[6] * y * z;
  res += coeff[7] * x * x;
  res += coeff[8] * y * y;
  res += coeff[9] * z * z;
  res += coeff[10] * x * x * y;
  res += coeff[11] * x * x * z;
  res += coeff[12] * y * y * x;
  res += coeff[13] * y * y * z;
  res += coeff[14] * z * z * x;
  res += coeff[15] * z * z * y;
  res += coeff[16] * x * x * x;
  res += coeff[17] * y * y * y;
  res += coeff[18] * z * z * z;
  res += coeff[19] * x * x * y * y;
  res += coeff[20] * x * x * z * z;
  res += coeff[21] * y * y * z * z;
  res += coeff[22] * x * x * x * y;
  res += coeff[23] * x * x * x * z;
  res += coeff[24] * y * y * y * x;
  res += coeff[25] * y * y * y * z;
  res += coeff[26] * z * z * z * x;
  res += coeff[27] * z * z * z * y;
  res += coeff[28] * x * x * x * x;
  res += coeff[29] * y * y * y * y;
  res += coeff[30] * z * z * z * z;
  return res;
}

double PolyApprox3D::val_at(const double x, const double y, const double z) {
  double res = 0.0;
  if(N == 2) {
    res = val_at_2nd(x - offset_x, y - offset_y, z - offset_z);
  } else if(N == 3) {
    res = val_at_3rd(x - offset_x, y - offset_y, z - offset_z);
  } else if(N == 4) {
    res = val_at_4th(x - offset_x, y - offset_y, z - offset_z);
  }
  return res;
}

void PolyApprox3D::grad_at_2nd(const double x, const double y, const double z,
                             double &dx, double &dy, double &dz) {
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;

  dx += coeff[1];
  dx += coeff[4] * y;
  dx += coeff[5] * z;
  dx += 2.0 * coeff[7] * x;

  dy += coeff[2];
  dy += coeff[4] * x;
  dy += coeff[6] * z;
  dy += 2.0 * coeff[8] * y;

  dz += coeff[3];
  dz += coeff[5] * x;
  dz += coeff[6] * y;
  dz += 2.0 * coeff[9] * z;
}

void PolyApprox3D::grad_at_3rd(const double x, const double y, const double z,
                             double &dx, double &dy, double &dz) {
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;

  dx += coeff[1];
  dx += coeff[4] * y;
  dx += coeff[5] * z;
  dx += 2.0 * coeff[7] * x;
  dx += 2.0 * coeff[10] * x * y;
  dx += 2.0 * coeff[11] * x * z;
  dx += coeff[12] * y * y;
  dx += coeff[14] * z * z;
  dx += 3.0 * coeff[16] * x * x;

  dy += coeff[2];
  dy += coeff[4] * x;
  dy += coeff[6] * z;
  dy += 2.0 * coeff[8] * y;
  dy += coeff[10] * x * x;
  dy += 2.0 * coeff[12] * y * x;
  dy += 2.0 * coeff[13] * y * z;
  dy += coeff[15] * z * z;
  dy += 3.0 * coeff[17] * y * y;

  dz += coeff[3];
  dz += coeff[5] * x;
  dz += coeff[6] * y;
  dz += 2.0 * coeff[9] * z;
  dz += coeff[11] * x * x;
  dz += coeff[13] * y * y;
  dz += 2.0 * coeff[14] * z * x;
  dz += 2.0 * coeff[15] * z * y;
  dz += 3.0 * coeff[18] * z * z;
}

void PolyApprox3D::grad_at_4th(const double x, const double y, const double z,
                             double &dx, double &dy, double &dz) {
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;

  dx += coeff[1];
  dx += coeff[4] * y;
  dx += coeff[5] * z;
  dx += 2.0 * coeff[7] * x;
  dx += 2.0 * coeff[10] * x * y;
  dx += 2.0 * coeff[11] * x * z;
  dx += coeff[12] * y * y;
  dx += coeff[14] * z * z;
  dx += 3.0 * coeff[16] * x * x;
  dx += 2.0 * coeff[19] * x * y * y;
  dx += 2.0 * coeff[20] * x * z * z;
  dx += 3.0 * coeff[22] * x * x * y;
  dx += 3.0 * coeff[23] * x * x * z;
  dx += coeff[24] * y * y * y;
  dx += coeff[26] * z * z * z;
  dx += 4.0 * coeff[28] * x * x * x;

  dy += coeff[2];
  dy += coeff[4] * x;
  dy += coeff[6] * z;
  dy += 2.0 * coeff[8] * y;
  dy += coeff[10] * x * x;
  dy += 2.0 * coeff[12] * y * x;
  dy += 2.0 * coeff[13] * y * z;
  dy += coeff[15] * z * z;
  dy += 3.0 * coeff[17] * y * y;
  dy += 2.0 * coeff[19] * x * x * y;
  dy += 2.0 * coeff[21] * y * z * z;
  dy += coeff[22] * x * x * x;
  dy += 3.0 * coeff[24] * y * y * x;
  dy += 3.0 * coeff[25] * y * y * z;
  dy += coeff[27] * z * z * z;
  dy += 4.0 * coeff[29] * y * y * y;

  dz += coeff[3];
  dz += coeff[5] * x;
  dz += coeff[6] * y;
  dz += 2.0 * coeff[9] * z;
  dz += coeff[11] * x * x;
  dz += coeff[13] * y * y;
  dz += 2.0 * coeff[14] * z * x;
  dz += 2.0 * coeff[15] * z * y;
  dz += 3.0 * coeff[18] * z * z;
  dz += 2.0 * coeff[20] * x * x * z;
  dz += 2.0 * coeff[21] * y * y * z;
  dz += coeff[23] * x * x * x;
  dz += coeff[25] * y * y * y;
  dz += 3.0 * coeff[26] * z * z * x;
  dz += 3.0 * coeff[27] * z * z * y;
  dz += 4.0 * coeff[30] * z * z * z;
}

void PolyApprox3D::grad_at(const double x, const double y, const double z,
                         double &dx, double &dy, double &dz) {
  if(N == 2) {
    grad_at_2nd(x - offset_x, y - offset_y, z - offset_z, dx, dy, dz);
  } else if(N == 3) {
    grad_at_3rd(x - offset_x, y - offset_y, z - offset_z, dx, dy, dz);
  } else if(N == 4) {
    grad_at_4th(x - offset_x, y - offset_y, z - offset_z, dx, dy, dz);
  }
}

void PolyApprox3D::hessian_at_2nd(const double x, const double y, const double z,
                                double &dx2, double &dy2, double &dz2,
                                double &dxy, double &dxz, double &dyz) {
  dx2 = 0.0;
  dy2 = 0.0;
  dz2 = 0.0;
  dxy = 0.0;
  dxz = 0.0;
  dyz = 0.0;

  dx2 += 2.0 * coeff[7];
  dy2 += 2.0 * coeff[8];
  dz2 += 2.0 * coeff[9];
  dxy += coeff[4];
  dxz += coeff[5];
  dyz += coeff[6];
}

void PolyApprox3D::hessian_at_3rd(const double x, const double y, const double z,
                                double &dx2, double &dy2, double &dz2,
                                double &dxy, double &dxz, double &dyz) {
  dx2 = 0.0;
  dy2 = 0.0;
  dz2 = 0.0;
  dxy = 0.0;
  dxz = 0.0;
  dyz = 0.0;

  dx2 += 2.0 * coeff[7];
  dx2 += 2.0 * coeff[10] * y;
  dx2 += 2.0 * coeff[11] * z;
  dx2 += 6.0 * coeff[16] * x;

  dy2 += 2.0 * coeff[8];
  dy2 += 2.0 * coeff[12] * x;
  dy2 += 2.0 * coeff[13] * z;
  dy2 += 6.0 * coeff[17] * y;

  dz2 += 2.0 * coeff[9];
  dz2 += 2.0 * coeff[14] * x;
  dz2 += 2.0 * coeff[15] * y;
  dz2 += 6.0 * coeff[18] * z;

  dxy += coeff[4];
  dxy += 2.0 * coeff[10] * x;
  dxy += 2.0 * coeff[12] * y;

  dxz += coeff[5];
  dxz += 2.0 * coeff[11] * x;
  dxz += 2.0 * coeff[14] * z;

  dyz += coeff[6];
  dyz += 2.0 * coeff[13] * y;
  dyz += 2.0 * coeff[15] * z;
}

void PolyApprox3D::hessian_at_4th(const double x, const double y, const double z,
                                double &dx2, double &dy2, double &dz2,
                                double &dxy, double &dxz, double &dyz) {
  dx2 = 0.0;
  dy2 = 0.0;
  dz2 = 0.0;
  dxy = 0.0;
  dxz = 0.0;
  dyz = 0.0;

  dx2 += 2.0 * coeff[7];
  dx2 += 2.0 * coeff[10] * y;
  dx2 += 2.0 * coeff[11] * z;
  dx2 += 6.0 * coeff[16] * x;
  dx2 += 2.0 * coeff[19] * y * y;
  dx2 += 2.0 * coeff[20] * z * z;
  dx2 += 6.0 * coeff[22] * x * y;
  dx2 += 6.0 * coeff[23] * x * z;
  dx2 += 12.0 * coeff[28] * x * x;

  dy2 += 2.0 * coeff[8];
  dy2 += 2.0 * coeff[12] * x;
  dy2 += 2.0 * coeff[13] * z;
  dy2 += 6.0 * coeff[17] * y;
  dy2 += 2.0 * coeff[19] * x * x;
  dy2 += 2.0 * coeff[21] * z * z;
  dy2 += 6.0 * coeff[24] * y * x;
  dy2 += 6.0 * coeff[25] * y * z;
  dy2 += 12.0 * coeff[29] * y * y;

  dz2 += 2.0 * coeff[9];
  dz2 += 2.0 * coeff[14] * x;
  dz2 += 2.0 * coeff[15] * y;
  dz2 += 6.0 * coeff[18] * z;
  dz2 += 2.0 * coeff[20] * x * x;
  dz2 += 2.0 * coeff[21] * y * y;
  dz2 += 6.0 * coeff[26] * z * x;
  dz2 += 6.0 * coeff[27] * z * y;
  dz2 += 12.0 * coeff[30] * z * z;

  dxy += coeff[4];
  dxy += 2.0 * coeff[10] * x;
  dxy += 2.0 * coeff[12] * y;
  dxy += 4.0 * coeff[19] * x * y;
  dxy += 3.0 * coeff[22] * x * x;
  dxy += 3.0 * coeff[24] * y * y;

  dxz += coeff[5];
  dxz += 2.0 * coeff[11] * x;
  dxz += 2.0 * coeff[14] * z;
  dxz += 4.0 * coeff[20] * x * z;
  dxz += 3.0 * coeff[23] * x * x;
  dxz += 3.0 * coeff[26] * z * z;

  dyz += coeff[6];
  dyz += 2.0 * coeff[13] * y;
  dyz += 2.0 * coeff[15] * z;
  dyz += 4.0 * coeff[21] * y * z;
  dyz += 3.0 * coeff[25] * y * y;
  dyz += 3.0 * coeff[27] * z * z;
}

void PolyApprox3D::hessian_at(const double x, const double y, const double z,
                            double &dx2, double &dy2, double &dz2,
                            double &dxy, double &dxz, double &dyz) {
  if(N == 2) {
    hessian_at_2nd(x - offset_x, y - offset_y, z - offset_z, dx2, dy2, dz2, dxy, dxz, dyz);
  } else if(N == 3) {
    hessian_at_3rd(x - offset_x, y - offset_y, z - offset_z, dx2, dy2, dz2, dxy, dxz, dyz);
  } else if(N == 4) {
    hessian_at_4th(x - offset_x, y - offset_y, z - offset_z, dx2, dy2, dz2, dxy, dxz, dyz);
  }
}

void num_pts_pos_neg(const vector<double> &s, int &pos, int &neg) {
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
    return 19;
  } else if(N == 4) {
    return 31;
  } else {
    return -1;
  }
}

int PolyApprox3D::num_pts() {
  if(N == 2) {
    return 10;
  } else if(N == 3) {
    return 50;
  } else if(N == 4) {
    return 31;
  } else {
    return 0;
  }
}

int PolyApprox3D::num_elem_stencil() {
  if(N == 2) {
    return 12;
  } else if(N == 3) {
    return 8;
  } else if(N == 4) {
    return 8;
  } else {
    return 0;
  }
}

double PolyApprox3D::get_coeff(int ind) {
  return coeff[ind];
}

void PolyApprox3D::get_offsets(double &x, double &y, double &z) {
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
  timer->endTimer("PolyApprox3D - get_stencils");
  return stencils;
}