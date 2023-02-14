#ifndef __INS_LS_REINIT_POLY_EVAL_CUDA_H
#define __INS_LS_REINIT_POLY_EVAL_CUDA_H

#include "dg_compiler_defs.h"

#include <vector>

#include "ls_reinit_poly.h"

class PolyEval {
public:
  __host__ PolyEval(std::vector<PolyApprox> &polys);
  __host__ ~PolyEval();

__device__ DG_FP val_at(const int cell_ind, const DG_FP x, const DG_FP y) {
  DG_FP res = 0.0;
  res += coeff[cell_ind * 10 + 0];
  res += coeff[cell_ind * 10 + 1] * x;
  res += coeff[cell_ind * 10 + 2] * y;
  res += coeff[cell_ind * 10 + 3] * x * x;
  res += coeff[cell_ind * 10 + 4] * x * y;
  res += coeff[cell_ind * 10 + 5] * y * y;
  res += coeff[cell_ind * 10 + 6] * x * x * x;
  res += coeff[cell_ind * 10 + 7] * x * x * y;
  res += coeff[cell_ind * 10 + 8] * x * y * y;
  res += coeff[cell_ind * 10 + 9] * y * y * y;
  return res;
}

__device__ void grad_at(const int cell_ind, const DG_FP x, const DG_FP y,
                                  DG_FP &dx, DG_FP &dy) {
  dx = 0.0;
  dy = 0.0;

  dx += coeff[cell_ind * 10 + 1];
  dx += 2.0 * coeff[cell_ind * 10 + 3] * x;
  dx += coeff[cell_ind * 10 + 4] * y;
  dx += 3.0 * coeff[cell_ind * 10 + 6] * x * x;
  dx += 2.0 * coeff[cell_ind * 10 + 7] * x * y;
  dx += coeff[cell_ind * 10 + 8] * y * y;

  dy += coeff[cell_ind * 10 + 2];
  dy += coeff[cell_ind * 10 + 4] * x;
  dy += 2.0 * coeff[cell_ind * 10 + 5] * y;
  dy += coeff[cell_ind * 10 + 7] * x * x;
  dy += 2.0 * coeff[cell_ind * 10 + 8] * x * y;
  dy += 3.0 * coeff[cell_ind * 10 + 9] * y * y;
}
__device__ void hessian_at(const int cell_ind, const DG_FP x, const DG_FP y,
                                     DG_FP &dx2, DG_FP &dxy, DG_FP &dy2) {
  dx2  = 2.0 * coeff[cell_ind * 10 + 3];
  dx2 += 6.0 * coeff[cell_ind * 10 + 6] * x;
  dx2 += 2.0 * coeff[cell_ind * 10 + 7] * y;

  dxy  = coeff[cell_ind * 10 + 4];
  dxy += 2.0 * coeff[cell_ind * 10 + 7] * x;
  dxy += 2.0 * coeff[cell_ind * 10 + 8] * y;

  dy2  = 2.0 * coeff[cell_ind * 10 + 5];
  dy2 += 2.0 * coeff[cell_ind * 10 + 8] * x;
  dy2 += 6.0 * coeff[cell_ind * 10 + 9] * y;
}

private:
  DG_FP *coeff;
};

#endif
