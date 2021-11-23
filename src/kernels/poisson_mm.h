inline void poisson_mm(const double *factor, const double *mm, double *op1) {
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_NP; n++) {
      int ind = m + n * DG_NP;
      op1[ind] += *factor * mm[ind];
    }
  }
}
