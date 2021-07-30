inline void poisson_apply_bc(const int *bedgeNum, const double *op,
                             const double *bc, double *rhs) {
  int exInd = *bedgeNum * 4;

  for(int m = 0; m < 6; m++) {
    int ind = m * 4;
    for(int n = 0; n < 4; n++) {
      rhs[m] += op[ind + n] * bc[exInd + n];
    }
  }
}
