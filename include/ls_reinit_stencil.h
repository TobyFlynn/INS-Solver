#ifndef __INS_LS_REINIT_STENCIL_H
#define __INS_LS_REINIT_STENCIL_H

#include "op_seq.h"

#include <vector>

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

class PolyApprox {
public:
  PolyApprox(const int order, const int cell_ind, 
             op_map edge_map, op_dat x_dat, op_dat y_dat, op_dat s_dat);

  double val_at(const double x, const double y);
  void grad_at(const double x, const double y, double &dx, double &dy);
  void hessian_at(const double x, const double y, double &dx2, 
                  double &dxy, double &dy2);

private:
  double offset_x, offset_y;
  std::vector<double> coeff;

  void get_offset(const int ind, op_dat x_dat, op_dat y_dat);
  void stencil_data(const std::vector<int> &stencil, op_dat x_dat, op_dat y_dat, op_dat s_dat, 
                    std::vector<double> &x, std::vector<double> &y, std::vector<double> &s);
  void stencil_ind(const int central_ind, op_map edge_map, std::vector<int> &stencil);
};

#endif