#ifndef __INS_LS_REINIT_POLY_H
#define __INS_LS_REINIT_POLY_H

#include "op_seq.h"

#include <vector>
#include <map>
#include <set>

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

class PolyApprox {
public:
  PolyApprox(const int cell_ind, std::set<int> stencil, const double *x_ptr,
             const double *y_ptr, const double *s_ptr);
  PolyApprox(std::vector<double> &c, double off_x, double off_y);

  double val_at(const double x, const double y);
  void grad_at(const double x, const double y, double &dx, double &dy);
  void hessian_at(const double x, const double y, double &dx2, 
                  double &dxy, double &dy2);
  double get_coeff(int ind);
  void get_offsets(double &x, double &y);

  static const int N = 2;
  static int num_coeff();
  static int num_elem_stencil();
  static std::map<int,std::set<int>> get_stencils(const std::set<int> &central_inds, op_map edge_map);
private:
  double offset_x, offset_y;
  std::vector<double> coeff;

  void get_offset(const int ind, const double *x_ptr, const double *y_ptr);
  void stencil_data(const std::set<int> &stencil, const double *x_ptr, const double *y_ptr, 
                    const double *s_ptr, std::vector<double> &x, std::vector<double> &y, 
                    std::vector<double> &s);
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