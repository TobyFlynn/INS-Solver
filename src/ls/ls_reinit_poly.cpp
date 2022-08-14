#include "op_seq.h"

#include "ls_reinit_poly.h"

#include <map>
#include <iostream>

#include "utils.h"

using namespace std;

bool vec_contains(const int val, const vector<int> &vec) {
  for(int i = 0; i < vec.size(); i++) {
    if(val == vec[i])
      return true;
  }
  return false;
}

void PolyApprox::get_offset(const int ind, const double *x_ptr, const double *y_ptr) {
  // offset_x = x_ptr[ind * DG_NP];
  // offset_y = y_ptr[ind * DG_NP];
  offset_x = 0.0;
  offset_y = 0.0;
}

struct Coord {
  double x;
  double y;
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
        if(xCmp && yCmp) {
          return false;
        } else if(xCmp) {
          return a.y < b.y;
        } else {
          return a.x < b.x;
        }
    }
};

void PolyApprox::stencil_data(const set<int> &stencil, const double *x_ptr, 
                              const double *y_ptr, const double *s_ptr, 
                              vector<double> &x, vector<double> &y, 
                              vector<double> &s) {
  map<Coord, Point, cmpCoords> pointMap;

  for(const auto &sten : stencil) {
    for(int n = 0; n < 6; n++) {
      int ind = sten * 6 + n;

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

  for(auto const &p : pointMap) {
    x.push_back(p.second.coord.x);
    y.push_back(p.second.coord.y);
    s.push_back(p.second.val / (double)p.second.count);
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

PolyApprox::PolyApprox(const int cell_ind, set<int> stencil, 
                       const double *x_ptr, const double *y_ptr, 
                       const double *s_ptr) {
  get_offset(cell_ind, x_ptr, y_ptr);
  
  vector<double> x_vec, y_vec, s_vec;
  stencil_data(stencil, x_ptr, y_ptr, s_ptr, x_vec, y_vec, s_vec);

  // Make sure equal number of points on each side of the line
  int pts_needed = num_coeff();
  int pts_aim = num_pts();
  int pos_pts, neg_pts;
  num_pts_pos_neg(s_vec, pos_pts, neg_pts);
  while((x_vec.size() > pts_needed && pos_pts != neg_pts) || x_vec.size() > pts_aim) {
    // Find point furthest from the interface to discard
    int ind_discard;
    if(pos_pts > neg_pts) {
      double max = s_vec[0];
      ind_discard = 0;
      for(int i = 1; i < x_vec.size(); i++) {
        if(s_vec[i] > max) {
          max = s_vec[i];
          ind_discard = i;
        }
      }
    } else {
      double min = s_vec[0];
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

PolyApprox::PolyApprox(std::vector<double> &c, double off_x, double off_y) {
  offset_x = off_x;
  offset_y = off_y;
  for(int i = 0; i < num_coeff(); i++) {
    coeff.push_back(c[i]);
  }
}

void PolyApprox::set_2nd_order_coeff(const vector<double> &x, const vector<double> &y, const vector<double> &s) {
  arma::mat A(x.size(), 6);
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,0) = 1.0;
    A(i,1) = x[i];
    A(i,2) = y[i];
    A(i,3) = x[i] * x[i];
    A(i,4) = x[i] * y[i];
    A(i,5) = y[i] * y[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  for(int i = 0; i < 6; i++) {
    coeff.push_back(ans(i));
  }
}

void PolyApprox::set_3rd_order_coeff(const vector<double> &x, const vector<double> &y, const vector<double> &s) {
  arma::mat A(x.size(), 10);
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,0) = 1.0;
    A(i,1) = x[i];
    A(i,2) = y[i];
    A(i,3) = x[i] * x[i];
    A(i,4) = x[i] * y[i];
    A(i,5) = y[i] * y[i];
    A(i,6) = x[i] * x[i] * x[i];
    A(i,7) = x[i] * x[i] * y[i];
    A(i,8) = x[i] * y[i] * y[i];
    A(i,9) = y[i] * y[i] * y[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  for(int i = 0; i < 10; i++) {
    coeff.push_back(ans(i));
  }
}

void PolyApprox::set_4th_order_coeff(const vector<double> &x, const vector<double> &y, const vector<double> &s) {
  arma::mat A(x.size(), 15);
  arma::vec b(x.size());
  for(int i = 0; i < x.size(); i++) {
    A(i,0)  = 1.0;
    A(i,1)  = x[i];
    A(i,2)  = y[i];
    A(i,3)  = x[i] * x[i];
    A(i,4)  = x[i] * y[i];
    A(i,5)  = y[i] * y[i];
    A(i,6)  = x[i] * x[i] * x[i];
    A(i,7)  = x[i] * x[i] * y[i];
    A(i,8)  = x[i] * y[i] * y[i];
    A(i,9)  = y[i] * y[i] * y[i];
    A(i,10) = x[i] * x[i] * x[i] * x[i];
    A(i,11) = x[i] * x[i] * x[i] * y[i];
    A(i,12) = x[i] * x[i] * y[i] * y[i];
    A(i,13) = x[i] * y[i] * y[i] * y[i];
    A(i,14) = y[i] * y[i] * y[i] * y[i];

    b(i) = s[i];
  }

  arma::vec ans = arma::solve(A, b);
  for(int i = 0; i < 15; i++) {
    coeff.push_back(ans(i));
  }
}

double PolyApprox::val_at_2nd(const double x, const double y) {
  double res = 0.0;
  res += coeff[0];
  res += coeff[1] * x;
  res += coeff[2] * y;
  res += coeff[3] * x * x;
  res += coeff[4] * x * y;
  res += coeff[5] * y * y;
  return res;
}

double PolyApprox::val_at_3rd(const double x, const double y) {
  double res = 0.0;
  res += coeff[0];
  res += coeff[1] * x;
  res += coeff[2] * y;
  res += coeff[3] * x * x;
  res += coeff[4] * x * y;
  res += coeff[5] * y * y;
  res += coeff[6] * x * x * x;
  res += coeff[7] * x * x * y;
  res += coeff[8] * x * y * y;
  res += coeff[9] * y * y * y;
  return res;
}

double PolyApprox::val_at_4th(const double x, const double y) {
  double res = 0.0;
  res += coeff[0];
  res += coeff[1] * x;
  res += coeff[2] * y;
  res += coeff[3] * x * x;
  res += coeff[4] * x * y;
  res += coeff[5] * y * y;
  res += coeff[6] * x * x * x;
  res += coeff[7] * x * x * y;
  res += coeff[8] * x * y * y;
  res += coeff[9] * y * y * y;
  res += coeff[10] * x * x * x * x;
  res += coeff[11] * x * x * x * y;
  res += coeff[12] * x * x * y * y;
  res += coeff[13] * x * y * y * y;
  res += coeff[14] * y * y * y * y;
  return res;
}

double PolyApprox::val_at(const double x, const double y) {
  double res = 0.0;
  if(N == 2) {
    res = val_at_2nd(x - offset_x, y - offset_y);
  } else if(N == 3) {
    res = val_at_3rd(x - offset_x, y - offset_y);
  } else if(N == 4) {
    res = val_at_4th(x - offset_x, y - offset_y);
  }
  return res;
}

void PolyApprox::grad_at_2nd(const double x, const double y, double &dx, double &dy) {
  dx = 0.0;
  dy = 0.0;
  
  dx += coeff[1];
  dx += 2.0 * coeff[3] * x;
  dx += coeff[4] * y;

  dy += coeff[2];
  dy += coeff[4] * x;
  dy += 2.0 * coeff[5] * y;
}

void PolyApprox::grad_at_3rd(const double x, const double y, double &dx, double &dy) {
  dx = 0.0;
  dy = 0.0;
  
  dx += coeff[1];
  dx += 2.0 * coeff[3] * x;
  dx += coeff[4] * y;
  dx += 3.0 * coeff[6] * x * x;
  dx += 2.0 * coeff[7] * x * y;
  dx += coeff[8] * y * y;

  dy += coeff[2];
  dy += coeff[4] * x;
  dy += 2.0 * coeff[5] * y;
  dy += coeff[7] * x * x;
  dy += 2.0 * coeff[8] * x * y;
  dy += 3.0 * coeff[9] * y * y;
}

void PolyApprox::grad_at_4th(const double x, const double y, double &dx, double &dy) {
  dx = 0.0;
  dy = 0.0;
  
  dx += coeff[1];
  dx += 2.0 * coeff[3] * x;
  dx += coeff[4] * y;
  dx += 3.0 * coeff[6] * x * x;
  dx += 2.0 * coeff[7] * x * y;
  dx += coeff[8] * y * y;
  dx += 4.0 * coeff[10] * x * x * x;
  dx += 3.0 * coeff[11] * x * x * y;
  dx += 2.0 * coeff[12] * x * y * y;
  dx += coeff[13] * y * y * y;

  dy += coeff[2];
  dy += coeff[4] * x;
  dy += 2.0 * coeff[5] * y;
  dy += coeff[7] * x * x;
  dy += 2.0 * coeff[8] * x * y;
  dy += 3.0 * coeff[9] * y * y;
  dy += coeff[11] * x * x * x;
  dy += 2.0 * coeff[12] * x * x * y;
  dy += 3.0 * coeff[13] * x * y * y;
  dy += 4.0 * coeff[14] * y * y * y;
}

void PolyApprox::grad_at(const double x, const double y, 
                         double &dx, double &dy) {
  if(N == 2) {
    grad_at_2nd(x - offset_x, y - offset_y, dx, dy);
  } else if(N == 3) {
    grad_at_3rd(x - offset_x, y - offset_y, dx, dy);
  } else if(N == 4) {
    grad_at_4th(x - offset_x, y - offset_y, dx, dy);
  }
}

void PolyApprox::hessian_at_2nd(const double x, const double y, double &dx2, double &dxy, double &dy2) {
  dx2 = 2.0 * coeff[3];
  dxy = coeff[4];
  dy2 = 2.0 * coeff[5];
}

void PolyApprox::hessian_at_3rd(const double x, const double y, double &dx2, double &dxy, double &dy2) {
  dx2  = 2.0 * coeff[3];
  dx2 += 6.0 * coeff[6] * x;
  dx2 += 2.0 * coeff[7] * y;

  dxy  = coeff[4];
  dxy += 2.0 * coeff[7] * x;
  dxy += 2.0 * coeff[8] * y;

  dy2  = 2.0 * coeff[5];
  dy2 += 2.0 * coeff[8] * x;
  dy2 += 6.0 * coeff[9] * y;
}

void PolyApprox::hessian_at_4th(const double x, const double y, double &dx2, double &dxy, double &dy2) {
  dx2  = 2.0 * coeff[3];
  dx2 += 6.0 * coeff[6] * x;
  dx2 += 2.0 * coeff[7] * y;
  dx2 += 12.0 * coeff[10] * x * x;
  dx2 += 6.0 * x * y;
  dx2 += 2.0 * y * y;

  dxy  = coeff[4];
  dxy += 2.0 * coeff[7] * x;
  dxy += 2.0 * coeff[8] * y;
  dxy += 3.0 * coeff[11] * x * x;
  dxy += 4.0 * coeff[12] * x * y;
  dxy += 3.0 * coeff[13] * y * y;

  dy2  = 2.0 * coeff[5];
  dy2 += 2.0 * coeff[8] * x;
  dy2 += 6.0 * coeff[9] * y;
  dy2 += 2.0 * coeff[12] * x * x;
  dy2 += 6.0 * coeff[13] * x * y;
  dy2 += 12.0 * coeff[14] * y * y;
}

void PolyApprox::hessian_at(const double x, const double y, 
                            double &dx2, double &dxy, double &dy2) {
  if(N == 2) {
    hessian_at_2nd(x - offset_x, y - offset_y, dx2, dxy, dy2);
  } else if(N == 3) {
    hessian_at_3rd(x - offset_x, y - offset_y, dx2, dxy, dy2);
  } else if(N == 4) {
    hessian_at_4th(x - offset_x, y - offset_y, dx2, dxy, dy2);
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
    return 25;
  } else if(N == 4) {
    return 25;
  } else {
    return 0;
  }
}

int PolyApprox::num_elem_stencil() {
  if(N == 2) {
    return 4;
  } else if(N == 3) {
    return 10;
  } else if(N == 4) {
    return 10;
  } else {
    return 0;
  }
}

double PolyApprox::get_coeff(int ind) {
  return coeff[ind];
}

void PolyApprox::get_offsets(double &x, double &y) {
  x = offset_x;
  y = offset_y;
}

struct stencil_query {
  int ind;
  set<int> central_inds;
};

map<int,set<int>> PolyApprox::get_stencils(const set<int> &central_inds, op_map edge_map) {
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
            if(stencil_it->second.find(edge_map->map[i + 1]) != stencil_it->second.end() 
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
            if(stencil_it->second.find(edge_map->map[i - 1]) != stencil_it->second.end() 
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
  
  return stencils;
}