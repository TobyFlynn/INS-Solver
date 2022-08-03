#include "ls_reinit_poly_eval_cuda.h"

PolyEval::PolyEval(std::map<int,PolyApprox> &polyMap) {
  auto firstPoly = polyMap.begin();
  const int numCoeffPerPoly = firstPoly->second.num_coeff();
  cudaMallocManaged(&coeff, polyMap.size() * numCoeffPerPoly * sizeof(double));

  int i = 0;
  for(auto it = polyMap.begin(); it != polyMap.end(); it++) {
    indMap.insert({it->first, i});
    for(int j = 0; j < numCoeffPerPoly; j++) {
      coeff[i * numCoeffPerPoly + j] = it->second.get_coeff(j);
    }
    i++;
  }
}

PolyEval::~PolyEval() {
  cudaFree(coeff);
}

int PolyEval::convert_ind(const int cell_ind) {
  return indMap.at(cell_ind);
}
/*
// For now just assume 3rd order poly
__device__ double PolyEval::val_at(const int cell_ind, const double x, const double y) {
  double res = 0.0;
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

__device__ void PolyEval::grad_at(const int cell_ind, const double x, const double y, 
                                  double &dx, double &dy) {
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
__device__ void PolyEval::hessian_at(const int cell_ind, const double x, const double y, 
                                     double &dx2, double &dxy, double &dy2) {
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
*/
