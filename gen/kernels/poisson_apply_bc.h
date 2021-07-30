inline void poisson_apply_bc(const int *bedgeNum, const double *op,
                             const double *bc, double *rhs) {
  int exInd = *bedgeNum * 3;

  for(int m = 0; m < 3; m++) {
    int ind = m * 3;
    for(int n = 0; n < 3; n++) {
      rhs[m] += op[ind + n] * bc[exInd + n];
    }
  }
}
