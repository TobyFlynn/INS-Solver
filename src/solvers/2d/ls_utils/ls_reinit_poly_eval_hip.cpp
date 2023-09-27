#include "ls_utils/2d/ls_reinit_poly_eval_device.h"

PolyEval::PolyEval(std::vector<PolyApprox> &polys) {
  throw std::runtime_error("PolyEval::PolyEval not implemented for HIP yet");
  /*const int numCoeffPerPoly = PolyApprox::num_coeff();
  cudaMallocManaged(&coeff, polys.size() * numCoeffPerPoly * sizeof(DG_FP));

  int i = 0;
  for(auto &poly : polys) {
    for(int j = 0; j < numCoeffPerPoly; j++) {
      coeff[i * numCoeffPerPoly + j] = poly.get_coeff(j);
    }
    i++;
  }*/
}

PolyEval::~PolyEval() {
  // cudaFree(coeff);
}
