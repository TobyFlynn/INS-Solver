inline void poisson_mm(const double *factor, const double *mm, double *op1) {
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m * DG_NP + n;
      int colInd = n * DG_NP + m;
      op1[ind] += *factor * mm[colInd];
    }
  }
}
