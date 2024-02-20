#include "ls_utils/2d/ls_reinit_poly.h"

#include <random>
#include <map>

#include "op2_utils.h"
#include "timing.h"
#include "dg_constants/dg_constants_2d.h"
#include "dg_utils.h"
#include "dg_global_constants/dg_global_constants_2d.h"

extern Timing *timer;
extern DGConstants *constants;

using namespace std;

struct Coord {
  DG_FP x;
  DG_FP y;
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
        if(xCmp && yCmp) {
          return false;
        } else if(xCmp) {
          return a.y < b.y;
        } else {
          return a.x < b.x;
        }
    }
};

void avg_stencil(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                 const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                 const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap);
void avg_stencil_nodes(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                 const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                 const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap);
void avg_stencil_if_contain_interface(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                 const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                 const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap);
void central_avg(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                 const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                 const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap);
void proj_onto_interface(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                         const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                         const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap);
void random_pts(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap);

bool vec_contains(const int val, const vector<int> &vec) {
  for(int i = 0; i < vec.size(); i++) {
    if(val == vec[i])
      return true;
  }
  return false;
}

void PolyApprox::calc_offset(const int ind, const DG_FP *x_ptr, const DG_FP *y_ptr) {
  offset_x = 0.0;
  offset_y = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    offset_x += x_ptr[ind * DG_NP + i];
    offset_y += y_ptr[ind * DG_NP + i];
  }
  offset_x /= (DG_FP)DG_NP;
  offset_y /= (DG_FP)DG_NP;
}

void PolyApprox::stencil_data(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                              const DG_FP *y_ptr, const DG_FP *s_ptr, 
                              const DG_FP *modal_ptr, vector<DG_FP> &x, 
                              vector<DG_FP> &y, vector<DG_FP> &s) {
  // Setup random number generator for later
  
  map<Coord, Point, cmpCoords> pointMap;
  avg_stencil(cell_ind, stencil, x_ptr, y_ptr, s_ptr, modal_ptr, offset_x, offset_y, pointMap);
  
  // avg_stencil_nodes(cell_ind, stencil, x_ptr, y_ptr, s_ptr, modal_ptr, offset_x, offset_y, pointMap);

  // central_avg(cell_ind, stencil, x_ptr, y_ptr, s_ptr, modal_ptr, offset_x, offset_y, pointMap);

  // proj_onto_interface(cell_ind, stencil, x_ptr, y_ptr, s_ptr, modal_ptr, offset_x, offset_y, pointMap);

  // random_pts(cell_ind, stencil, x_ptr, y_ptr, s_ptr, modal_ptr, offset_x, offset_y, pointMap);

  for(auto const &p : pointMap) {
    x.push_back(p.second.coord.x);
    y.push_back(p.second.coord.y);
    s.push_back(p.second.val / (DG_FP)p.second.count);
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

PolyApprox::PolyApprox(const int cell_ind, set<int> stencil,
                       const DG_FP *x_ptr, const DG_FP *y_ptr,
                       const DG_FP *s_ptr, const DG_FP *modal_ptr) {
  calc_offset(cell_ind, x_ptr, y_ptr);

  vector<DG_FP> x_vec, y_vec, s_vec;
  stencil_data(cell_ind, stencil, x_ptr, y_ptr, s_ptr, modal_ptr, x_vec, y_vec, s_vec);

  // Make sure equal number of points on each side of the line
  int pts_needed = num_coeff();
  int pts_aim = num_pts();
  int pos_pts, neg_pts;
  num_pts_pos_neg(s_vec, pos_pts, neg_pts);
  while((x_vec.size() > pts_needed && pos_pts != neg_pts) || x_vec.size() > pts_aim) {
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
    s_vec.erase(s_vec.begin() + ind_discard);

    num_pts_pos_neg(s_vec, pos_pts, neg_pts);
  }

  if(N == 2) {
    set_2nd_order_coeff(x_vec, y_vec, s_vec);
  } else if(N == 3) {
    set_3rd_order_coeff(x_vec, y_vec, s_vec);
  } else if(N == 4) {
    set_4th_order_coeff(x_vec, y_vec, s_vec);
  }
}

PolyApprox::PolyApprox(std::vector<DG_FP> &c, DG_FP off_x, DG_FP off_y) {
  offset_x = off_x;
  offset_y = off_y;
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(c[i]);
  }
}

#define C_POLY_IND 0
#define X_POLY_IND 1
#define Y_POLY_IND 2
#define X2_POLY_IND 3
#define XY_POLY_IND 4
#define Y2_POLY_IND 5
#define X3_POLY_IND 6
#define X2Y_POLY_IND 7
#define Y2X_POLY_IND 8
#define Y3_POLY_IND 9
#define X4_POLY_IND 10
#define X3Y_POLY_IND 11
#define X2Y2_POLY_IND 12
#define Y3X_POLY_IND 13
#define Y4_POLY_IND 14

void PolyApprox::set_2nd_order_coeff(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &s) {
  arma::mat A(x.size(), 6);
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,C_POLY_IND)  = 1.0;
    A(i,X_POLY_IND)  = x[i];
    A(i,Y_POLY_IND)  = y[i];
    A(i,X2_POLY_IND) = x[i] * x[i];
    A(i,XY_POLY_IND) = x[i] * y[i];
    A(i,Y2_POLY_IND) = y[i] * y[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  // arma::vec ans = arma::solve(arma::inv(A.t() * A) * A.t(), b);
  for(int i = 0; i < 6; i++) {
    coeff.push_back(ans(i));
  }
}

void PolyApprox::set_3rd_order_coeff(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &s) {
  arma::mat A(x.size(), 10);
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,C_POLY_IND)   = 1.0;
    A(i,X_POLY_IND)   = x[i];
    A(i,Y_POLY_IND)   = y[i];
    A(i,X2_POLY_IND)  = x[i] * x[i];
    A(i,XY_POLY_IND)  = x[i] * y[i];
    A(i,Y2_POLY_IND)  = y[i] * y[i];
    A(i,X3_POLY_IND)  = x[i] * x[i] * x[i];
    A(i,X2Y_POLY_IND) = x[i] * x[i] * y[i];
    A(i,Y2X_POLY_IND) = x[i] * y[i] * y[i];
    A(i,Y3_POLY_IND)  = y[i] * y[i] * y[i];

    b(i) = s[i];
  }

  // arma::vec ans = arma::solve(A, b);
  arma::vec ans = arma::inv(A.t() * A) * A.t() * b;
  for(int i = 0; i < 10; i++) {
    coeff.push_back(ans(i));
  }
}

void PolyApprox::set_4th_order_coeff(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &s) {
  arma::mat A(x.size(), 15);
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,C_POLY_IND)    = 1.0;
    A(i,X_POLY_IND)    = x[i];
    A(i,Y_POLY_IND)    = y[i];
    A(i,X2_POLY_IND)   = x[i] * x[i];
    A(i,XY_POLY_IND)   = x[i] * y[i];
    A(i,Y2_POLY_IND)   = y[i] * y[i];
    A(i,X3_POLY_IND)   = x[i] * x[i] * x[i];
    A(i,X2Y_POLY_IND)  = x[i] * x[i] * y[i];
    A(i,Y2X_POLY_IND)  = x[i] * y[i] * y[i];
    A(i,Y3_POLY_IND)   = y[i] * y[i] * y[i];
    A(i,X4_POLY_IND)   = x[i] * x[i] * x[i] * x[i];
    A(i,X3Y_POLY_IND)  = x[i] * x[i] * x[i] * y[i];
    A(i,X2Y2_POLY_IND) = x[i] * x[i] * y[i] * y[i];
    A(i,Y3X_POLY_IND)  = x[i] * y[i] * y[i] * y[i];
    A(i,Y4_POLY_IND)   = y[i] * y[i] * y[i] * y[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  // arma::vec ans = arma::inv(A.t() * A) * A.t() * b;
  for(int i = 0; i < 15; i++) {
    coeff.push_back(ans(i));
  }
}

DG_FP PolyApprox::val_at_2nd(const DG_FP x, const DG_FP y) {
  DG_FP res = 0.0;
  res += coeff[C_POLY_IND];
  res += coeff[X_POLY_IND] * x;
  res += coeff[Y_POLY_IND] * y;
  res += coeff[X2_POLY_IND] * x * x;
  res += coeff[XY_POLY_IND] * x * y;
  res += coeff[Y2_POLY_IND] * y * y;
  return res;
}

DG_FP PolyApprox::val_at_3rd(const DG_FP x, const DG_FP y) {
  DG_FP res = 0.0;
  res += coeff[C_POLY_IND];
  res += coeff[X_POLY_IND] * x;
  res += coeff[Y_POLY_IND] * y;
  res += coeff[X2_POLY_IND] * x * x;
  res += coeff[XY_POLY_IND] * x * y;
  res += coeff[Y2_POLY_IND] * y * y;
  res += coeff[X3_POLY_IND] * x * x * x;
  res += coeff[X2Y_POLY_IND] * x * x * y;
  res += coeff[Y2X_POLY_IND] * x * y * y;
  res += coeff[Y3_POLY_IND] * y * y * y;
  return res;
}

DG_FP PolyApprox::val_at_4th(const DG_FP x, const DG_FP y) {
  DG_FP res = 0.0;
  res += coeff[C_POLY_IND];
  res += coeff[X_POLY_IND] * x;
  res += coeff[Y_POLY_IND] * y;
  res += coeff[X2_POLY_IND] * x * x;
  res += coeff[XY_POLY_IND] * x * y;
  res += coeff[Y2_POLY_IND] * y * y;
  res += coeff[X3_POLY_IND] * x * x * x;
  res += coeff[X2Y_POLY_IND] * x * x * y;
  res += coeff[Y2X_POLY_IND] * x * y * y;
  res += coeff[Y3_POLY_IND] * y * y * y;
  res += coeff[X4_POLY_IND] * x * x * x * x;
  res += coeff[X3Y_POLY_IND] * x * x * x * y;
  res += coeff[X2Y2_POLY_IND] * x * x * y * y;
  res += coeff[Y3X_POLY_IND] * x * y * y * y;
  res += coeff[Y4_POLY_IND] * y * y * y * y;
  return res;
}

DG_FP PolyApprox::val_at(const DG_FP x, const DG_FP y) {
  DG_FP res = 0.0;
  if(N == 2) {
    res = val_at_2nd(x, y);
  } else if(N == 3) {
    res = val_at_3rd(x, y);
  } else if(N == 4) {
    res = val_at_4th(x, y);
  }
  return res;
}

void PolyApprox::grad_at_2nd(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy) {
  dx = 0.0;
  dy = 0.0;

  dx += coeff[X_POLY_IND];
  dx += 2.0 * coeff[X2_POLY_IND] * x;
  dx += coeff[XY_POLY_IND] * y;

  dy += coeff[Y_POLY_IND];
  dy += coeff[XY_POLY_IND] * x;
  dy += 2.0 * coeff[Y2_POLY_IND] * y;
}

void PolyApprox::grad_at_3rd(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy) {
  dx = 0.0;
  dy = 0.0;

  dx += coeff[X_POLY_IND];
  dx += 2.0 * coeff[X2_POLY_IND] * x;
  dx += coeff[XY_POLY_IND] * y;
  dx += 3.0 * coeff[X3_POLY_IND] * x * x;
  dx += 2.0 * coeff[X2Y_POLY_IND] * x * y;
  dx += coeff[Y2X_POLY_IND] * y * y;

  dy += coeff[Y_POLY_IND];
  dy += coeff[XY_POLY_IND] * x;
  dy += 2.0 * coeff[Y2_POLY_IND] * y;
  dy += coeff[X2Y_POLY_IND] * x * x;
  dy += 2.0 * coeff[Y2X_POLY_IND] * x * y;
  dy += 3.0 * coeff[Y3_POLY_IND] * y * y;
}

void PolyApprox::grad_at_4th(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy) {
  dx = 0.0;
  dy = 0.0;

  dx += coeff[X_POLY_IND];
  dx += 2.0 * coeff[X2_POLY_IND] * x;
  dx += coeff[XY_POLY_IND] * y;
  dx += 3.0 * coeff[X3_POLY_IND] * x * x;
  dx += 2.0 * coeff[X2Y_POLY_IND] * x * y;
  dx += coeff[Y2X_POLY_IND] * y * y;
  dx += 4.0 * coeff[X4_POLY_IND] * x * x * x;
  dx += 3.0 * coeff[X3Y_POLY_IND] * x * x * y;
  dx += 2.0 * coeff[X2Y2_POLY_IND] * x * y * y;
  dx += coeff[Y3X_POLY_IND] * y * y * y;

  dy += coeff[Y_POLY_IND];
  dy += coeff[XY_POLY_IND] * x;
  dy += 2.0 * coeff[Y2_POLY_IND] * y;
  dy += coeff[X2Y_POLY_IND] * x * x;
  dy += 2.0 * coeff[Y2X_POLY_IND] * x * y;
  dy += 3.0 * coeff[Y3_POLY_IND] * y * y;
  dy += coeff[X3Y_POLY_IND] * x * x * x;
  dy += 2.0 * coeff[X2Y2_POLY_IND] * x * x * y;
  dy += 3.0 * coeff[Y3X_POLY_IND] * x * y * y;
  dy += 4.0 * coeff[Y4_POLY_IND] * y * y * y;
}

void PolyApprox::grad_at(const DG_FP x, const DG_FP y,
                         DG_FP &dx, DG_FP &dy) {
  if(N == 2) {
    grad_at_2nd(x, y, dx, dy);
  } else if(N == 3) {
    grad_at_3rd(x, y, dx, dy);
  } else if(N == 4) {
    grad_at_4th(x, y, dx, dy);
  }
}

void PolyApprox::hessian_at_2nd(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  dx2 = 2.0 * coeff[X2_POLY_IND];
  dxy = coeff[XY_POLY_IND];
  dy2 = 2.0 * coeff[Y2_POLY_IND];
}

void PolyApprox::hessian_at_3rd(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  dx2  = 2.0 * coeff[X2_POLY_IND];
  dx2 += 6.0 * coeff[X3_POLY_IND] * x;
  dx2 += 2.0 * coeff[X2Y_POLY_IND] * y;

  dxy  = coeff[XY_POLY_IND];
  dxy += 2.0 * coeff[X2Y_POLY_IND] * x;
  dxy += 2.0 * coeff[Y2X_POLY_IND] * y;

  dy2  = 2.0 * coeff[Y2_POLY_IND];
  dy2 += 2.0 * coeff[Y2X_POLY_IND] * x;
  dy2 += 6.0 * coeff[Y3_POLY_IND] * y;
}

void PolyApprox::hessian_at_4th(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  dx2  = 2.0 * coeff[X2_POLY_IND];
  dx2 += 6.0 * coeff[X3_POLY_IND] * x;
  dx2 += 2.0 * coeff[X2Y_POLY_IND] * y;
  dx2 += 12.0 * coeff[X4_POLY_IND] * x * x;
  dx2 += 6.0 * coeff[X3Y_POLY_IND] * x * y;
  dx2 += 2.0 * coeff[X2Y2_POLY_IND] * y * y;

  dxy  = coeff[XY_POLY_IND];
  dxy += 2.0 * coeff[X2Y_POLY_IND] * x;
  dxy += 2.0 * coeff[Y2X_POLY_IND] * y;
  dxy += 3.0 * coeff[X3Y_POLY_IND] * x * x;
  dxy += 4.0 * coeff[X2Y2_POLY_IND] * x * y;
  dxy += 3.0 * coeff[Y3X_POLY_IND] * y * y;

  dy2  = 2.0 * coeff[Y2_POLY_IND];
  dy2 += 2.0 * coeff[Y2X_POLY_IND] * x;
  dy2 += 6.0 * coeff[Y3_POLY_IND] * y;
  dy2 += 2.0 * coeff[X2Y2_POLY_IND] * x * x;
  dy2 += 6.0 * coeff[Y3X_POLY_IND] * x * y;
  dy2 += 12.0 * coeff[Y4_POLY_IND] * y * y;
}

void PolyApprox::hessian_at(const DG_FP x, const DG_FP y,
                            DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  if(N == 2) {
    hessian_at_2nd(x, y, dx2, dxy, dy2);
  } else if(N == 3) {
    hessian_at_3rd(x, y, dx2, dxy, dy2);
  } else if(N == 4) {
    hessian_at_4th(x, y, dx2, dxy, dy2);
  }
}

int PolyApprox::num_coeff() {
  if(N == 2) {
    return 6;
  } else if(N == 3) {
    return 10;
  } else if(N == 4) {
    return 15;
  } else {
    return -1;
  }
}

int PolyApprox::num_pts() {
  if(N == 2) {
    return 12;
  } else if(N == 3) {
    return 26;
  } else if(N == 4) {
    return 25;
  } else {
    return 0;
  }
}

int PolyApprox::num_elem_stencil() {
  if(N == 2) {
    return 13;
  } else if(N == 3) {
    return 13;
  } else if(N == 4) {
    return 13;
  } else {
    return 0;
  }
}

DG_FP PolyApprox::get_coeff(int ind) {
  return coeff[ind];
}

void PolyApprox::get_offsets(DG_FP &x, DG_FP &y) {
  x = offset_x;
  y = offset_y;
}

struct stencil_query {
  int ind;
  set<int> central_inds;
};

map<int,set<int>> PolyApprox::get_stencils(const set<int> &central_inds, op_map edge_map, const DG_FP *x_ptr, const DG_FP *y_ptr) {
  return single_layer_stencils(central_inds, edge_map, x_ptr, y_ptr);
  
  timer->startTimer("PolyApprox - get_stencils");
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
  timer->endTimer("PolyApprox - get_stencils");

  return stencils;
}

bool share_coords(const DG_FP *x_ptr, const DG_FP *y_ptr, const std::vector<Coord> &nodes) {
  for(int i = 0; i < DG_NP; i++) {
    for(int n = 0; n < 3; n++) {
      bool xCmp = abs(x_ptr[i] - nodes[n].x) < 1e-8;
      bool yCmp = abs(y_ptr[i] - nodes[n].y) < 1e-8;
      if(xCmp && yCmp) return true;
    }
  }
  return false;
}

map<int,set<int>> PolyApprox::single_layer_stencils(const set<int> &central_inds, op_map edge_map, const DG_FP *x_ptr, const DG_FP *y_ptr) {
  timer->startTimer("PolyApprox - get_stencils");
  map<int,set<int>> stencils;
  map<int,stencil_query> queryInds;
  const int num_elements = num_elem_stencil();
  map<int,std::vector<Coord>> central_inds_nodes;

  const int fmask_node_ind_0 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int fmask_node_ind_1 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + DG_NPF - 1];
  const int fmask_node_ind_2 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + 2 * DG_NPF - 1];

  for(const auto &ind : central_inds) {
    set<int> st;
    st.insert(ind);
    stencils.insert({ind, st});
    stencil_query sq;
    sq.ind = ind;
    sq.central_inds.insert(ind);
    queryInds.insert({ind, sq});

    std::vector<Coord> nodes;
    Coord node0, node1, node2;
    node0.x = x_ptr[ind * DG_NP + fmask_node_ind_0];
    node0.y = y_ptr[ind * DG_NP + fmask_node_ind_0];
    nodes.push_back(node0);
    node1.x = x_ptr[ind * DG_NP + fmask_node_ind_1];
    node1.y = y_ptr[ind * DG_NP + fmask_node_ind_1];
    nodes.push_back(node1);
    node2.x = x_ptr[ind * DG_NP + fmask_node_ind_2];
    node2.y = y_ptr[ind * DG_NP + fmask_node_ind_2];
    nodes.push_back(node2);
    central_inds_nodes.insert({ind, nodes});
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
              // Check if we share a node with the central ind
              auto node_coords = central_inds_nodes.at(ind);
              if(share_coords(x_ptr + edge_map->map[i + 1] * DG_NP, y_ptr + edge_map->map[i + 1] * DG_NP, node_coords)) {
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
          }
        } else {
          // For each central ind associated with this query ind
          for(const auto &ind : it->second.central_inds) {
            auto stencil_it = stencils.find(ind);
            // Check if the other cell in this edge is already in the stencil for this central ind
            if(stencil_it->second.find(edge_map->map[i - 1]) == stencil_it->second.end()
               && stencil_it->second.size() < num_elements) {
              auto node_coords = central_inds_nodes.at(ind);
              if(share_coords(x_ptr + edge_map->map[i - 1] * DG_NP, y_ptr + edge_map->map[i - 1] * DG_NP, node_coords)) {
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
    }

    queryInds = newQueryInds;
  }
  timer->endTimer("PolyApprox - get_stencils");

  return stencils;
}

/*
 * Different ways of getting the stencil data
 */

DG_FP pol_in_tri_sign(const DG_FP x, const DG_FP y, const DG_FP v0x, const DG_FP v0y, const DG_FP v1x, const DG_FP v1y) {
  return (x - v1x) * (v0y - v1y) - (v0x - v1x) * (y - v1y);
}

bool pol_in_tri(const DG_FP ptX, const DG_FP ptY, const DG_FP *nodeX, const DG_FP *nodeY) {
  bool has_neg, has_pos;

  DG_FP d1 = pol_in_tri_sign(ptX, ptY, nodeX[0], nodeY[0], nodeX[1], nodeY[1]);
  DG_FP d2 = pol_in_tri_sign(ptX, ptY, nodeX[1], nodeY[1], nodeX[2], nodeY[2]);
  DG_FP d3 = pol_in_tri_sign(ptX, ptY, nodeX[2], nodeY[2], nodeX[0], nodeY[0]);

  has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

  return !(has_neg && has_pos);
}

void pol_rs2xy(DG_FP &sampleX, DG_FP &sampleY, const DG_FP *nodeX, const DG_FP *nodeY) {
  DG_FP r_ = sampleX;
  DG_FP s_ = sampleY;

  sampleX = 0.5 * (nodeX[1] * (1.0 + r_) + nodeX[2] * (1.0 + s_) - nodeX[0] * (s_ + r_));
  sampleY = 0.5 * (nodeY[1] * (1.0 + r_) + nodeY[2] * (1.0 + s_) - nodeY[0] * (s_ + r_));
}

bool pol_simplified_newton(DG_FP &pt_r, DG_FP &pt_s, const DG_FP *modal,
                           const DG_FP tol, const DG_FP *tetra_r, const DG_FP *tetra_s) {
  bool converged = false;
  for(int step = 0; step < 50; step++) {
    DG_FP surf = DGUtils::val_at_pt_2d(pt_r, pt_s, modal, DG_ORDER);
    DG_FP dsdr, dsds;
    DGUtils::grad_at_pt_2d(pt_r, pt_s, modal, DG_ORDER, dsdr, dsds);

    DG_FP sqrnorm = dsdr * dsdr + dsds * dsds;
    if(sqrnorm > 0.0) {
      dsdr *= surf / sqrnorm;
      dsds *= surf / sqrnorm;
    }

    pt_r -= dsdr;
    pt_s -= dsds;

    if(!pol_in_tri(pt_r, pt_s, tetra_r, tetra_s)) {
      break;
    }

    // Check convergence
    if(dsdr * dsdr + dsds * dsds < tol) {
      converged = true;
      break;
    }
  }
  return converged;
}

void avg_stencil(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                 const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                 const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap) {
  for(const auto &sten : stencil) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = sten * DG_NP + n;

      Coord coord;
      coord.x = x_ptr[ind] - offset_x;
      coord.y = y_ptr[ind] - offset_y;
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
}

void avg_stencil_nodes(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                 const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                 const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap) {
  const int fmask_node_ind_0 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int fmask_node_ind_1 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + DG_NPF - 1];
  const int fmask_node_ind_2 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + 2 * DG_NPF - 1];
  int inds[] = {fmask_node_ind_0, fmask_node_ind_1, fmask_node_ind_2};
  for(const auto &sten : stencil) {

    for(int n = 0; n < 3; n++) {
      int ind = sten * DG_NP + inds[n];

      Coord coord;
      coord.x = x_ptr[ind] - offset_x;
      coord.y = y_ptr[ind] - offset_y;
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
}

void avg_stencil_if_contain_interface(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                 const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                 const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap) {
  for(const auto &sten : stencil) {
    bool contain_interface = false;
    bool pos = s_ptr[sten * DG_NP] > 0.0;
    for(int i = 1; i < DG_NP; i++) {
      if(pos != s_ptr[sten * DG_NP + i] > 0.0)
        contain_interface = true;
    }
    if(!contain_interface) continue;
    for(int n = 0; n < DG_NP; n++) {
      int ind = sten * DG_NP + n;

      Coord coord;
      coord.x = x_ptr[ind] - offset_x;
      coord.y = y_ptr[ind] - offset_y;
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
}

void central_avg(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                 const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                 const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap) {
  for(int n = 0; n < DG_NP; n++) {
    int ind = cell_ind * DG_NP + n;
    Coord coord;
    coord.x = x_ptr[ind] - offset_x;
    coord.y = y_ptr[ind] - offset_y;
    Point point;
    auto res = pointMap.insert(make_pair(coord, point));
    res.first->second.coord = coord;
    res.first->second.val   = s_ptr[ind];
    res.first->second.count = 1;
  }

  for(const auto &sten : stencil) {
    if(sten == cell_ind) continue;
    for(int n = 0; n < DG_NP; n++) {
      int ind = sten * DG_NP + n;

      Coord coord;
      coord.x = x_ptr[ind] - offset_x;
      coord.y = y_ptr[ind] - offset_y;
      Point point;
      auto res = pointMap.find(coord);
      if(res == pointMap.end()) continue;

      res->second.val += s_ptr[ind];
      res->second.count++;
    }
  }
}

void proj_onto_interface(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                         const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                         const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1.0, 1.0);

  const int fmask_node_ind_0 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int fmask_node_ind_1 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + DG_NPF - 1];
  const int fmask_node_ind_2 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + 2 * DG_NPF - 1];
  const DG_FP *tmp_r_ptr = constants->get_mat_ptr(DGConstants::R) + (DG_ORDER - 1) * DG_NP;
  const DG_FP *tmp_s_ptr = constants->get_mat_ptr(DGConstants::S) + (DG_ORDER - 1) * DG_NP;
  const DG_FP tetra_r[3] = {tmp_r_ptr[fmask_node_ind_0], tmp_r_ptr[fmask_node_ind_1], tmp_r_ptr[fmask_node_ind_2]};
  const DG_FP tetra_s[3] = {tmp_s_ptr[fmask_node_ind_0], tmp_s_ptr[fmask_node_ind_1], tmp_s_ptr[fmask_node_ind_2]};

  for(const auto &sten : stencil) {
    bool contain_interface = false;
    bool s_pos = s_ptr[sten * DG_NP] > 0.0;
    for(int i = 0; i < DG_NP; i++) {
      if(s_ptr[sten * DG_NP + i] > 0.0 != s_pos) {
        contain_interface = true;
        break;
      }
    }
    if(!contain_interface || sten != cell_ind) continue;

    DG_FP ref_r_ptr[DG_NP], ref_s_ptr[DG_NP];
    for(int i = 0; i < DG_NP; i++) {
      ref_r_ptr[i] = tmp_r_ptr[i] * 0.75;
      ref_s_ptr[i] = tmp_s_ptr[i] * 0.75;
    }

    const DG_FP X0 = x_ptr[sten * DG_NP + fmask_node_ind_0]; const DG_FP Y0 = y_ptr[sten * DG_NP + fmask_node_ind_0];
    const DG_FP X1 = x_ptr[sten * DG_NP + fmask_node_ind_1]; const DG_FP Y1 = y_ptr[sten * DG_NP + fmask_node_ind_1];
    const DG_FP X2 = x_ptr[sten * DG_NP + fmask_node_ind_2]; const DG_FP Y2 = y_ptr[sten * DG_NP + fmask_node_ind_2];

    const DG_FP nodeX[] = {X0, X1, X2};
    const DG_FP nodeY[] = {Y0, Y1, Y2};

    #pragma omp parallel for
    for(int p = 0; p < DG_NP; p++) {
      // bool converged = false;
      bool converged = pol_simplified_newton(ref_r_ptr[p], ref_s_ptr[p], &modal_ptr[sten * DG_NP], 1e-8, tetra_r, tetra_s);

      if(!converged && sten == cell_ind) {
        int counter = 0;
        while(!converged && counter < 10) {
          ref_r_ptr[p] = dis(gen);
          ref_s_ptr[p] = dis(gen);
          while(!pol_in_tri(ref_r_ptr[p], ref_s_ptr[p], tetra_r, tetra_s)) {
            ref_r_ptr[p] = dis(gen);
            ref_s_ptr[p] = dis(gen);
          }
          converged = pol_simplified_newton(ref_r_ptr[p], ref_s_ptr[p], &modal_ptr[sten * DG_NP], 1e-8, tetra_r, tetra_s);
          counter++;
        }
      }

      // Check if point has converged
      // && in_tetra(ref_r_ptr[p], ref_s_ptr[p], ref_t_ptr[p], tetra_r, tetra_s, tetra_t)
      if(converged) {
        DG_FP x = ref_r_ptr[p];
        DG_FP y = ref_s_ptr[p];
        pol_rs2xy(x, y, nodeX, nodeY);
        const DG_FP tmp_val = DGUtils::val_at_pt_2d(ref_r_ptr[p], ref_s_ptr[p], &modal_ptr[sten * DG_NP], DG_ORDER);
        #pragma omp critical
        {
          Coord coord;
          coord.x = x - offset_x;
          coord.y = y - offset_y;
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
}

void random_pts(const int cell_ind, const set<int> &stencil, const DG_FP *x_ptr,
                const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr, 
                const DG_FP offset_x, const DG_FP offset_y, map<Coord, Point, cmpCoords> &pointMap) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1.0, 1.0);

  const int fmask_node_ind_0 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  const int fmask_node_ind_1 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + DG_NPF - 1];
  const int fmask_node_ind_2 = FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + 2 * DG_NPF - 1];
  const DG_FP X0 = x_ptr[cell_ind * DG_NP + fmask_node_ind_0]; const DG_FP Y0 = y_ptr[cell_ind * DG_NP + fmask_node_ind_0];
  const DG_FP X1 = x_ptr[cell_ind * DG_NP + fmask_node_ind_1]; const DG_FP Y1 = y_ptr[cell_ind * DG_NP + fmask_node_ind_1];
  const DG_FP X2 = x_ptr[cell_ind * DG_NP + fmask_node_ind_2]; const DG_FP Y2 = y_ptr[cell_ind * DG_NP + fmask_node_ind_2];
  const DG_FP nodeX[] = {X0, X1, X2};
  const DG_FP nodeY[] = {Y0, Y1, Y2};

  for(int p = 0; p < DG_NP; p++) {
    DG_FP sampleX = dis(gen);
    DG_FP sampleY = dis(gen);
    DG_FP surf = DGUtils::val_at_pt_2d(sampleX, sampleY, &modal_ptr[cell_ind * DG_NP], DG_ORDER);
    pol_rs2xy(sampleX, sampleY, nodeX, nodeY);
    while(!pol_in_tri(sampleX, sampleY, nodeX, nodeY)) {
      sampleX = dis(gen);
      sampleY = dis(gen);
      surf = DGUtils::val_at_pt_2d(sampleX, sampleY, &modal_ptr[cell_ind * DG_NP], DG_ORDER);
      pol_rs2xy(sampleX, sampleY, nodeX, nodeY);
    }

    Coord coord;
    coord.x = sampleX - offset_x;
    coord.y = sampleY - offset_y;
    Point point;
    auto res = pointMap.insert(make_pair(coord, point));
    if(res.second) {
      // Point was inserted
      res.first->second.coord = coord;
      res.first->second.val   = surf;
      res.first->second.count = 1;
    }
  }
}