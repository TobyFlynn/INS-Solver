#include "shared_save_functions.h"

#include <memory>
#include <algorithm>
#include <map>
#include <cmath>

using namespace std;

struct Point {
  double x;
  double y;
  double u;
  double v;
  double pr;
  double vort;
  double s;
  vector<int> cells;
  vector<int> pointNum;
  int counter;
};

struct cmpCoords {
    bool operator()(const pair<double,double>& a, const pair<double,double>& b) const {
        bool xCmp = abs(a.first - b.first) < 1e-8;
        bool yCmp = abs(a.second - b.second) < 1e-8;
        if(xCmp && yCmp) {
          return false;
        }
        return a < b;
    }
};

void get_data_vectors_order_4(vector<double> &x_v, vector<double> &y_v,
                              vector<double> &u_v, vector<double> &v_v,
                              vector<double> &pr_v, vector<double> &vort_v,
                              vector<double> &s_v, vector<cgsize_t> &cells,
                              double *Ux, double *Uy, double *pr, double *vort,
                              double *x, double *y, double *s, int numCells) {
  // Maps points to sub elements that they are part of.
  // Each line is 6 long (as 6 is the max number of sub elements within an original element that a point can be part of)
  // -1 is just padding to get each line to 6
  int cellMask[15][6] = {
    {0, -1, -1, -1, -1, -1}, // Point 0 is part of sub element 0
    {0, 1, 2, -1, -1, -1},    // 1
    {2, 3, 4, -1, -1, -1},    // 2
    {4, 5, 6, -1, -1, -1},    // 3
    {6, -1, -1, -1, -1, -1}, // End of first point row
    {0, 1, 7, -1, -1, -1},    // 5
    {1, 2, 3, 9, 8, 7},       // 6
    {3, 4, 5, 11, 10, 9},     // 7
    {6, 5, 11, -1, -1, -1}, // End of second point row
    {7, 8, 12, -1, -1, -1},   // 9
    {8, 9, 10, 14, 13, 12},   // 10
    {11, 10, 14, -1, -1, -1}, // End of third point row
    {12, 13, 15, -1, -1, -1}, // 12
    {14, 13, 15, -1, -1, -1}, // End of fourth point row
    {15, -1, -1, -1, -1, -1}  // 14
  };

  int pointNum[15][6] = {
    {0, -1, -1, -1, -1, -1}, // Point 0 is the first point of sub element 0
    {1, 0, 0, -1, -1, -1},  // 1
    {1, 0, 0, -1, -1, -1},  // 2
    {1, 0, 0, -1, -1, -1},  // 3
    {1, -1, -1, -1, -1, -1}, // End of first point row
    {2, 2, 0, -1, -1, -1},  // 5
    {1, 2, 2, 0, 0, 1},     // 6
    {1, 2, 2, 0, 0, 1},     // 7
    {2, 1, 1, -1, -1, -1}, // End of second point row
    {2, 2, 0, -1, -1, -1},  // 9
    {1, 2, 2, 0, 0, 1},     // 10
    {2, 1, 1, -1, -1, -1}, // End of third point row
    {2, 2, 0, -1, -1, -1},  // 12
    {2, 1, 1, -1, -1, -1}, // End of fourth point row
    {2, -1, -1, -1, -1, -1} // 14
  };

  // 16 sub elements per original element
  map<pair<double,double>,unique_ptr<Point>, cmpCoords> pointMap;
  for(int c = 0; c < numCells; c++) {
    int ind = c * 15;
    for(int p = 0; p < 15; p++) {
      pair<double,double> coords = make_pair(x[ind + p], y[ind + p]);
      unique_ptr<Point> point = make_unique<Point>();
      auto res = pointMap.insert(make_pair(coords, move(point)));
      if(res.second) {
        res.first->second->x    = x[ind + p];
        res.first->second->y    = y[ind + p];
        res.first->second->u    = Ux[ind + p];
        res.first->second->v    = Uy[ind + p];
        res.first->second->pr   = pr[ind + p];
        res.first->second->vort = vort[ind + p];
        res.first->second->s    = s[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 16 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter = 1;
      } else {
        res.first->second->u    += Ux[ind + p];
        res.first->second->v    += Uy[ind + p];
        res.first->second->pr   += pr[ind + p];
        res.first->second->vort += vort[ind + p];
        res.first->second->s    += s[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 16 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter++;
      }
    }
  }

  int index = 0;

  for(auto const &p : pointMap) {
    x_v.push_back(p.second->x);
    y_v.push_back(p.second->y);
    u_v.push_back(p.second->u / p.second->counter);
    v_v.push_back(p.second->v / p.second->counter);
    pr_v.push_back(p.second->pr / p.second->counter);
    vort_v.push_back(p.second->vort / p.second->counter);
    s_v.push_back(p.second->s / p.second->counter);
    for(int i = 0; i < p.second->cells.size(); i++) {
      cells[p.second->cells[i] * 3 + p.second->pointNum[i]] = index + 1;
    }
    index++;
  }
}

void get_data_vectors_order_3(vector<double> &x_v, vector<double> &y_v,
                              vector<double> &u_v, vector<double> &v_v,
                              vector<double> &pr_v, vector<double> &vort_v,
                              vector<double> &s_v, vector<cgsize_t> &cells,
                              double *Ux, double *Uy, double *pr, double *vort,
                              double *x, double *y, double *s, int numCells) {
  // Maps points to sub elements that they are part of.
  // Each line is 6 long (as 6 is the max number of sub elements within an original element that a point can be part of)
  // -1 is just padding to get each line to 6
  int cellMask[10][6] = {
    {0, -1, -1, -1, -1, -1},  // Point 0 is part of sub element 0
    {0, 1, 2, -1, -1, -1},    // 1
    {2, 3, 4, -1, -1, -1},    // 2
    {4, -1, -1, -1, -1, -1},  // End of first point row
    {0, 1, 5, -1, -1, -1},    // 4
    {1, 2, 3, 7, 6, 5},       // 5
    {4, 3, 7, -1, -1, -1},    // End of second point row
    {5, 6, 8, -1, -1, -1},    // 7
    {7, 6, 8, -1, -1, -1},    // End of third point row
    {8, -1, -1, -1, -1, -1}   // 9
  };

  int pointNum[10][6] = {
    {0, -1, -1, -1, -1, -1},  // Point 0 is the first point of sub element 0
    {1, 0, 0, -1, -1, -1},    // 1
    {1, 0, 0, -1, -1, -1},    // 2
    {1, -1, -1, -1, -1, -1},  // End of first point row
    {2, 2, 0, -1, -1, -1},    // 4
    {1, 2, 2, 0, 0, 1},       // 5
    {2, 1, 1, -1, -1, -1},    // End of second point row
    {2, 2, 0, -1, -1, -1},    // 7
    {2, 1, 1, -1, -1, -1},    // End of third point row
    {2, -1, -1, -1, -1, -1}   // 9
  };

  // 9 sub elements per original element
  map<pair<double,double>,unique_ptr<Point>, cmpCoords> pointMap;
  for(int c = 0; c < numCells; c++) {
    int ind = c * 3;
    for(int p = 0; p < 3; p++) {
      pair<double,double> coords = make_pair(x[ind + p], y[ind + p]);
      unique_ptr<Point> point = make_unique<Point>();
      auto res = pointMap.insert(make_pair(coords, move(point)));
      if(res.second) {
        res.first->second->x    = x[ind + p];
        res.first->second->y    = y[ind + p];
        res.first->second->u    = Ux[ind + p];
        res.first->second->v    = Uy[ind + p];
        res.first->second->pr   = pr[ind + p];
        res.first->second->vort = vort[ind + p];
        res.first->second->s    = s[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 9 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter = 1;
      } else {
        res.first->second->u    += Ux[ind + p];
        res.first->second->v    += Uy[ind + p];
        res.first->second->pr   += pr[ind + p];
        res.first->second->vort += vort[ind + p];
        res.first->second->s    += s[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 9 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter++;
      }
    }
  }

  int index = 0;

  for(auto const &p : pointMap) {
    x_v.push_back(p.second->x);
    y_v.push_back(p.second->y);
    u_v.push_back(p.second->u / p.second->counter);
    v_v.push_back(p.second->v / p.second->counter);
    pr_v.push_back(p.second->pr / p.second->counter);
    vort_v.push_back(p.second->vort / p.second->counter);
    s_v.push_back(p.second->s / p.second->counter);
    for(int i = 0; i < p.second->cells.size(); i++) {
      cells[p.second->cells[i] * 3 + p.second->pointNum[i]] = index + 1;
    }
    index++;
  }
}

void get_data_vectors_order_2(vector<double> &x_v, vector<double> &y_v,
                              vector<double> &u_v, vector<double> &v_v,
                              vector<double> &pr_v, vector<double> &vort_v,
                              vector<double> &s_v, vector<cgsize_t> &cells,
                              double *Ux, double *Uy, double *pr, double *vort,
                              double *x, double *y, double *s, int numCells) {
  // Maps points to sub elements that they are part of.
  // Each line is 6 long (as 6 is the max number of sub elements within an original element that a point can be part of)
  // -1 is just padding to get each line to 6
  int cellMask[6][6] = {
    {0, -1, -1, -1, -1, -1},  // Point 0 is part of sub element 0
    {0, 1, 2, -1, -1, -1},    // 1
    {2, -1, -1, -1, -1, -1},  // End of first point row
    {0, 1, 3, -1, -1, -1},    // 3
    {2, 1, 3, -1, -1, -1},    // End of second point row
    {3, -1, -1, -1, -1, -1}   // 5
  };

  int pointNum[6][6] = {
    {0, -1, -1, -1, -1, -1},  // Point 0 is the first point of sub element 0
    {1, 0, 0, -1, -1, -1},    // 1
    {1, -1, -1, -1, -1, -1},  // End of first point row
    {2, 2, 0, -1, -1, -1},    // 3
    {2, 1, 1, -1, -1, -1},    // End of second point row
    {2, -1, -1, -1, -1, -1}   // 5
  };

  // 4 sub elements per original element
  map<pair<double,double>,unique_ptr<Point>, cmpCoords> pointMap;
  for(int c = 0; c < numCells; c++) {
    int ind = c * 3;
    for(int p = 0; p < 3; p++) {
      pair<double,double> coords = make_pair(x[ind + p], y[ind + p]);
      unique_ptr<Point> point = make_unique<Point>();
      auto res = pointMap.insert(make_pair(coords, move(point)));
      if(res.second) {
        res.first->second->x    = x[ind + p];
        res.first->second->y    = y[ind + p];
        res.first->second->u    = Ux[ind + p];
        res.first->second->v    = Uy[ind + p];
        res.first->second->pr   = pr[ind + p];
        res.first->second->vort = vort[ind + p];
        res.first->second->s    = s[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 4 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter = 1;
      } else {
        res.first->second->u    += Ux[ind + p];
        res.first->second->v    += Uy[ind + p];
        res.first->second->pr   += pr[ind + p];
        res.first->second->vort += vort[ind + p];
        res.first->second->s    += s[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 4 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter++;
      }
    }
  }

  int index = 0;

  for(auto const &p : pointMap) {
    x_v.push_back(p.second->x);
    y_v.push_back(p.second->y);
    u_v.push_back(p.second->u / p.second->counter);
    v_v.push_back(p.second->v / p.second->counter);
    pr_v.push_back(p.second->pr / p.second->counter);
    vort_v.push_back(p.second->vort / p.second->counter);
    s_v.push_back(p.second->s / p.second->counter);
    for(int i = 0; i < p.second->cells.size(); i++) {
      cells[p.second->cells[i] * 3 + p.second->pointNum[i]] = index + 1;
    }
    index++;
  }
}

void get_data_vectors_order_1(vector<double> &x_v, vector<double> &y_v,
                              vector<double> &u_v, vector<double> &v_v,
                              vector<double> &pr_v, vector<double> &vort_v,
                              vector<double> &s_v, vector<cgsize_t> &cells,
                              double *Ux, double *Uy, double *pr, double *vort,
                              double *x, double *y, double *s, int numCells) {
  // Maps points to sub elements that they are part of.
  // Each line is 6 long (as 6 is the max number of sub elements within an original element that a point can be part of)
  // -1 is just padding to get each line to 6
  int cellMask[3][6] = {
    {0, -1, -1, -1, -1, -1},  // Point 0 is part of sub element 0
    {0, -1, -1, -1, -1, -1},  // 1
    {0, -1, -1, -1, -1, -1}
  };

  int pointNum[3][6] = {
    {0, -1, -1, -1, -1, -1},  // Point 0 is the first point of sub element 0
    {1, -1, -1, -1, -1, -1},  // 1
    {2, -1, -1, -1, -1, -1}
  };

  // 4 sub elements per original element
  map<pair<double,double>,unique_ptr<Point>, cmpCoords> pointMap;
  for(int c = 0; c < numCells; c++) {
    int ind = c * 3;
    for(int p = 0; p < 3; p++) {
      pair<double,double> coords = make_pair(x[ind + p], y[ind + p]);
      unique_ptr<Point> point = make_unique<Point>();
      auto res = pointMap.insert(make_pair(coords, move(point)));
      if(res.second) {
        res.first->second->x    = x[ind + p];
        res.first->second->y    = y[ind + p];
        res.first->second->u    = Ux[ind + p];
        res.first->second->v    = Uy[ind + p];
        res.first->second->pr   = pr[ind + p];
        res.first->second->vort = vort[ind + p];
        res.first->second->s    = s[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 1 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter = 1;
      } else {
        res.first->second->u    += Ux[ind + p];
        res.first->second->v    += Uy[ind + p];
        res.first->second->pr   += pr[ind + p];
        res.first->second->vort += vort[ind + p];
        res.first->second->s    += s[ind + p];
        for(int m = 0; m < 6; m++) {
          if(cellMask[p][m] >= 0) {
            res.first->second->cells.push_back(c * 1 + cellMask[p][m]);
            res.first->second->pointNum.push_back(pointNum[p][m]);
          } else {
            break;
          }
        }
        res.first->second->counter++;
      }
    }
  }

  int index = 0;

  for(auto const &p : pointMap) {
    x_v.push_back(p.second->x);
    y_v.push_back(p.second->y);
    u_v.push_back(p.second->u / p.second->counter);
    v_v.push_back(p.second->v / p.second->counter);
    pr_v.push_back(p.second->pr / p.second->counter);
    vort_v.push_back(p.second->vort / p.second->counter);
    s_v.push_back(p.second->s / p.second->counter);
    for(int i = 0; i < p.second->cells.size(); i++) {
      cells[p.second->cells[i] * 3 + p.second->pointNum[i]] = index + 1;
    }
    index++;
  }
}
