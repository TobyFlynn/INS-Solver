inline void poisson_apply_bc(const int *bedgeNum, const double *op,
                             const double *bc, double *rhs) {
  int exInd = *bedgeNum * DG_GF_NP;

  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_GF_NP; n++) {
      int ind = m + n * DG_NP;
      rhs[m] += op[ind] * bc[exInd + n];
    }
  }
}
