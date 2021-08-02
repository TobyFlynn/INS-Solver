inline void poisson_apply_bc(const int *bedgeNum, const double *op,
                             const double *bc, double *rhs) {
  int exInd = *bedgeNum * 6;

  for(int m = 0; m < 10; m++) {
    int ind = m * 6;
    for(int n = 0; n < 6; n++) {
      rhs[m] += op[ind + n] * bc[exInd + n];
    }
  }
}
