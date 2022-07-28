#include "op_seq.h"

#include "ls_reinit_stencil.h"

#include <map>
#include <iostream>

using namespace std;

bool vec_contains(const int val, const vector<int> &vec) {
  for(int i = 0; i < vec.size(); i++) {
    if(val == vec[i])
      return true;
  }
  return false;
}

void PolyApprox::get_offset(const int ind, op_dat x_dat, op_dat y_dat) {
  const double *x_ptr = (double *)x_dat->data;
  const double *y_ptr = (double *)y_dat->data;
  // offset_x = x_ptr[ind * DG_NP];
  // offset_y = y_ptr[ind * DG_NP];
  offset_x = 0.0;
  offset_y = 0.0;
}

// Get vector of indices that make up stencil around central_ind
void PolyApprox::stencil_ind(const int central_ind, op_map edge_map, vector<int> &stencil) {
  vector<int> sweep1;
  const int numEdges = edge_map->from->size;
  for(int i = 0; i < numEdges * 2; i++) {
    if(edge_map->map[i] == central_ind) {
      if(i % 2 == 0) {
        sweep1.push_back(edge_map->map[i + 1]);
      } else {
        sweep1.push_back(edge_map->map[i - 1]);
      }

      if(sweep1.size() == 3)
        break;
    }
  }

  // vector<int> sweep2;
  // for(int i = 0; i < numEdges * 2; i++) {
  //   if(vec_contains(edge_map->map[i], sweep1)) {
  //     if(i % 2 == 0) {
  //       if(edge_map->map[i + 1] != central_ind)
  //         sweep2.push_back(edge_map->map[i + 1]);
  //     } else {
  //       if(edge_map->map[i - 1] != central_ind)
  //         sweep2.push_back(edge_map->map[i - 1]);
  //     }

  //     if(sweep2.size() == 6)
  //       break;
  //   }
  // }

  stencil.push_back(central_ind);

  for(int i = 0; i < sweep1.size(); i++) {
    stencil.push_back(sweep1[i]);
  }

  // for(int i = 0; i < sweep2.size(); i++) {
  //   stencil.push_back(sweep2[i]);
  // }
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

void PolyApprox::stencil_data(const vector<int> &stencil, op_dat x_dat, op_dat y_dat, 
                              op_dat s_dat, vector<double> &x, vector<double> &y, 
                              vector<double> &s) {
  map<Coord, Point, cmpCoords> pointMap;
  // TODO surround in OP2 MPI exchange thing if necessary
  const double *x_ptr = (double *)x_dat->data;
  const double *y_ptr = (double *)y_dat->data;
  const double *s_ptr = (double *)s_dat->data;

  for(int i = 0; i < stencil.size(); i++) {
    for(int n = 0; n < 6; n++) {
      int ind = stencil[i] * 6 + n;

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

void num_pts_pos_neg(const vector<double> s, int &pos, int &neg) {
  pos = 0;
  neg = 0;
  for(int i = 0; i < s.size(); i++) {
    if(s[i] > 0.0) pos++;
    if(s[i] < 0.0) neg++;
  }
}

PolyApprox::PolyApprox(const int order, const int cell_ind, op_map edge_map, 
                       op_dat x_dat, op_dat y_dat, op_dat s_dat) {
  get_offset(cell_ind, x_dat, y_dat);
  
  vector<int> stencil;
  stencil_ind(cell_ind, edge_map, stencil);
  vector<double> x_vec, y_vec, s_vec;
  stencil_data(stencil, x_dat, y_dat, s_dat, x_vec, y_vec, s_vec);

  // Make sure equal number of points on each side of the line
  int pos_pts, neg_pts;
  num_pts_pos_neg(s_vec, pos_pts, neg_pts);
  while(x_vec.size() > 6 && pos_pts != neg_pts) {
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

  // Start with second order poly
  arma::mat A(x_vec.size(), 6);
  arma::vec b(x_vec.size());
  for(int i = 0; i < x_vec.size(); i++) {
    A(i,0) = 1.0;
    A(i,1) = x_vec[i];
    A(i,2) = y_vec[i];
    A(i,3) = x_vec[i] * x_vec[i];
    A(i,4) = x_vec[i] * y_vec[i];
    A(i,5) = y_vec[i] * y_vec[i];

    b(i) = s_vec[i];
  }

  arma::vec x = arma::solve(A, b);
  for(int i = 0; i < 6; i++) {
    coeff.push_back(x(i));
  }
}

double PolyApprox::val_at(const double x, const double y) {
  double res = 0.0;
  res += coeff[0];
  res += coeff[1] * x;
  res += coeff[2] * y;
  res += coeff[3] * x * x;
  res += coeff[4] * x * y;
  res += coeff[5] * y * y;
  return res;
}

void PolyApprox::grad_at(const double x, const double y, 
                         double &dx, double &dy) {
  dx = 0.0;
  dy = 0.0;
  
  dx += coeff[1];
  dx += 2.0 * coeff[3] * x;
  dx += coeff[4] * y;

  dy += coeff[2];
  dy += coeff[4] * x;
  dy += 2.0 * coeff[5] * y;
}

void PolyApprox::hessian_at(const double x, const double y, 
                            double &dx2, double &dxy, double &dy2) {
  dx2 = 2.0 * coeff[3];
  dxy = coeff[4];
  dy2 = 2.0 * coeff[5];
}