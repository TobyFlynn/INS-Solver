#ifndef __INS_LS_REINIT_POLY_EVAL_CUDA_H
#define __INS_LS_REINIT_POLY_EVAL_CUDA_H

#include <map>

#include "ls_reinit_poly.h"

class PolyEval {
public:
  __host__ PolyEval(std::map<int,PolyApprox> &polyMap);
  __host__ ~PolyEval();

  __host__ int convert_ind(const int cell_ind);

  __device__ double val_at(const int cell_ind, const double x, const double y);
  __device__ void grad_at(const int cell_ind, const double x, const double y, 
                          double &dx, double &dy);
  __device__ void hessian_at(const int cell_ind, const double x, const double y, 
                             double &dx2, double &dxy, double &dy2);

private:
  double *coeff;
  std::map<int,int> indMap;
};

#endif