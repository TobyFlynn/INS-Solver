inline void poisson_mf2_apply_bc(const int *bedgeNum, const double *op,
                                 const double *bc, double *rhs) {
  int exInd = 0;
  if(*bedgeNum == 1) exInd = 7;
  else if(*bedgeNum == 2) exInd = 14;
  
  for(int m = 0; m < 15; m++) {
    int ind = m * 7;
    double val = 0.0;
    for(int n = 0; n < 7; n++) {
      val += op[ind + n] * bc[exInd + n];
    }
    rhs[m] += val;
  }
}
