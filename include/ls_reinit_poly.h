#ifndef __INS_LS_REINIT_POLY_H
#define __INS_LS_REINIT_POLY_H

#include "op_seq.h"

#include <vector>

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

class PolyApprox {
public:
  PolyApprox(const int order, const int cell_ind, 
             op_map edge_map, const double *x_ptr,
             const double *y_ptr, const double *s_ptr);

  double val_at(const double x, const double y);
  void grad_at(const double x, const double y, double &dx, double &dy);
  void hessian_at(const double x, const double y, double &dx2, 
                  double &dxy, double &dy2);
  int num_coeff();
  double get_coeff(int ind);
private:
  double offset_x, offset_y;
  std::vector<double> coeff;
  int N;

  void get_offset(const int ind, const double *x_ptr, const double *y_ptr);
  void stencil_data(const std::vector<int> &stencil, const double *x_ptr, const double *y_ptr, 
                    const double *s_ptr, std::vector<double> &x, std::vector<double> &y, 
                    std::vector<double> &s);
  void stencil_ind(const int central_ind, const int num_sweeps, op_map edge_map, 
                   std::vector<int> &stencil);
  int num_pts();

  void set_2nd_order_coeff(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s);
  void set_3rd_order_coeff(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s);
  void set_4th_order_coeff(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &s);

  double val_at_2nd(const double x, const double y);
  double val_at_3rd(const double x, const double y);
  double val_at_4th(const double x, const double y);

  void grad_at_2nd(const double x, const double y, double &dx, double &dy);
  void grad_at_3rd(const double x, const double y, double &dx, double &dy);
  void grad_at_4th(const double x, const double y, double &dx, double &dy);

  void hessian_at_2nd(const double x, const double y, double &dx2, double &dxy, double &dy2);
  void hessian_at_3rd(const double x, const double y, double &dx2, double &dxy, double &dy2);
  void hessian_at_4th(const double x, const double y, double &dx2, double &dxy, double &dy2);
};

#endif