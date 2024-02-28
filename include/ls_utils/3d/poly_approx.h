#ifndef __INS_POLY_APPROX_3D_H
#define __INS_POLY_APPROX_3D_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include <vector>
#include <map>
#include <set>

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

class PolyApprox3D {
public:
  PolyApprox3D(const int cell_ind, std::set<int> stencil, const DG_FP *x_ptr,
               const DG_FP *y_ptr, const DG_FP *z_ptr, const DG_FP *s_ptr, const DG_FP h);
  PolyApprox3D(std::vector<DG_FP> &c, DG_FP off_x, DG_FP off_y, DG_FP off_z);

  DG_FP val_at(const DG_FP x, const DG_FP y, const DG_FP z);
  void grad_at(const DG_FP x, const DG_FP y, const DG_FP z,
               DG_FP &dx, DG_FP &dy, DG_FP &dz);
  void hessian_at(const DG_FP x, const DG_FP y, const DG_FP z,
                  DG_FP &dx2, DG_FP &dy2, DG_FP &dz2,
                  DG_FP &dxy, DG_FP &dxz, DG_FP &dyz);
  DG_FP get_coeff(int ind);
  void get_offsets(DG_FP &x, DG_FP &y, DG_FP &z);

  static const int N = DG_ORDER;
  static int num_coeff();
  static int num_elem_stencil();
  static std::map<int,std::set<int>> get_stencils(const std::set<int> &central_inds, op_map edge_map);
private:
  DG_FP offset_x, offset_y, offset_z;
  std::vector<DG_FP> coeff;

  void get_offset(const int ind, const DG_FP *x_ptr, const DG_FP *y_ptr, const DG_FP *z_ptr);
  void stencil_data(const int cell_ind, const std::set<int> &stencil, const DG_FP *x_ptr, const DG_FP *y_ptr,
                    const DG_FP *z_ptr, const DG_FP *s_ptr, std::vector<DG_FP> &x,
                    std::vector<DG_FP> &y, std::vector<DG_FP> &z, std::vector<DG_FP> &s);
  int num_pts();

  void fit_poly(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &z, const vector<DG_FP> &s, const DG_FP h);
  arma::mat get_vandermonde(const vector<DG_FP> &x, const vector<DG_FP> &y, const vector<DG_FP> &z);
  arma::mat get_2nd_order_vandermonde(const std::vector<DG_FP> &x, const std::vector<DG_FP> &y, const std::vector<DG_FP> &z);
  arma::mat get_3rd_order_vandermonde(const std::vector<DG_FP> &x, const std::vector<DG_FP> &y, const std::vector<DG_FP> &z);
  arma::mat get_4th_order_vandermonde(const std::vector<DG_FP> &x, const std::vector<DG_FP> &y, const std::vector<DG_FP> &z);

  DG_FP val_at_2nd(const DG_FP x, const DG_FP y, const DG_FP z);
  DG_FP val_at_3rd(const DG_FP x, const DG_FP y, const DG_FP z);
  DG_FP val_at_4th(const DG_FP x, const DG_FP y, const DG_FP z);

  void grad_at_2nd(const DG_FP x, const DG_FP y, const DG_FP z,
                   DG_FP &dx, DG_FP &dy, DG_FP &dz);
  void grad_at_3rd(const DG_FP x, const DG_FP y, const DG_FP z,
                   DG_FP &dx, DG_FP &dy, DG_FP &dz);
  void grad_at_4th(const DG_FP x, const DG_FP y, const DG_FP z,
                   DG_FP &dx, DG_FP &dy, DG_FP &dz);

  void hessian_at_2nd(const DG_FP x, const DG_FP y, const DG_FP z,
                      DG_FP &dx2, DG_FP &dy2, DG_FP &dz2,
                      DG_FP &dxy, DG_FP &dxz, DG_FP &dyz);
  void hessian_at_3rd(const DG_FP x, const DG_FP y, const DG_FP z,
                      DG_FP &dx2, DG_FP &dy2, DG_FP &dz2,
                      DG_FP &dxy, DG_FP &dxz, DG_FP &dyz);
  void hessian_at_4th(const DG_FP x, const DG_FP y, const DG_FP z,
                      DG_FP &dx2, DG_FP &dy2, DG_FP &dz2,
                      DG_FP &dxy, DG_FP &dxz, DG_FP &dyz);
};

#endif
