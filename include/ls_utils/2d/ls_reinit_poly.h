#ifndef __INS_LS_REINIT_POLY_H
#define __INS_LS_REINIT_POLY_H

#include "dg_compiler_defs.h"

#include "op_seq.h"

#include <vector>
#include <map>
#include <set>

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

class PolyApprox {
public:
  PolyApprox(const int cell_ind, std::set<int> stencil, const DG_FP *x_ptr,
             const DG_FP *y_ptr, const DG_FP *s_ptr, const DG_FP *modal_ptr);
  PolyApprox(std::vector<DG_FP> &c, DG_FP off_x, DG_FP off_y);

  DG_FP val_at(const DG_FP x, const DG_FP y);
  void grad_at(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy);
  void hessian_at(const DG_FP x, const DG_FP y, DG_FP &dx2,
                  DG_FP &dxy, DG_FP &dy2);
  DG_FP get_coeff(int ind);
  void get_offsets(DG_FP &x, DG_FP &y);

  static const int N = DG_ORDER;
  static int num_coeff();
  static int num_elem_stencil();
  static std::map<int,std::set<int>> get_stencils(const std::set<int> &central_inds, op_map edge_map, const DG_FP *x_ptr, const DG_FP *y_ptr);
private:
  DG_FP offset_x, offset_y;
  std::vector<DG_FP> coeff;

  static std::map<int,std::set<int>> single_layer_stencils(const std::set<int> &central_inds, op_map edge_map, const DG_FP *x_ptr, const DG_FP *y_ptr);

  void calc_offset(const int ind, const DG_FP *x_ptr, const DG_FP *y_ptr);
  void stencil_data(const int cell_ind, const std::set<int> &stencil, const DG_FP *x_ptr, const DG_FP *y_ptr,
                    const DG_FP *s_ptr, const DG_FP *modal_ptr, std::vector<DG_FP> &x, std::vector<DG_FP> &y,
                    std::vector<DG_FP> &s);
  int num_pts();

  void set_2nd_order_coeff(const std::vector<DG_FP> &x, const std::vector<DG_FP> &y, const std::vector<DG_FP> &s);
  void set_3rd_order_coeff(const std::vector<DG_FP> &x, const std::vector<DG_FP> &y, const std::vector<DG_FP> &s);
  void set_4th_order_coeff(const std::vector<DG_FP> &x, const std::vector<DG_FP> &y, const std::vector<DG_FP> &s);

  DG_FP val_at_2nd(const DG_FP x, const DG_FP y);
  DG_FP val_at_3rd(const DG_FP x, const DG_FP y);
  DG_FP val_at_4th(const DG_FP x, const DG_FP y);

  void grad_at_2nd(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy);
  void grad_at_3rd(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy);
  void grad_at_4th(const DG_FP x, const DG_FP y, DG_FP &dx, DG_FP &dy);

  void hessian_at_2nd(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2);
  void hessian_at_3rd(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2);
  void hessian_at_4th(const DG_FP x, const DG_FP y, DG_FP &dx2, DG_FP &dxy, DG_FP &dy2);
};

#endif
