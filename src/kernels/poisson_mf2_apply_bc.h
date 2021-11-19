inline void poisson_mf2_apply_bc(const int *bedgeNum, const double *op,
                                 const double *bc, double *rhs) {
  int exInd = *bedgeNum * DG_GF_NP;

  for(int m = 0; m < DG_NP; m++) {
    int ind = m * DG_GF_NP;
    double val = 0.0;
    for(int n = 0; n < DG_GF_NP; n++) {
      val += op[ind + n] * bc[exInd + n];
    }
    rhs[m] += val;
  }
}
