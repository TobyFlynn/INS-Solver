inline void poisson_apply_bc(const int *bedgeNum, const double *op,
                             const double *bc, double *rhs) {
  int exInd = *bedgeNum * 7;

  for(int m = 0; m < 15; m++) {
    int ind = m * 7;
    for(int n = 0; n < 7; n++) {
      rhs[m] += op[ind + n] * bc[exInd + n];
    }
  }
}
