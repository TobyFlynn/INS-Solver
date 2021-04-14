inline void poisson_test_bc(const int *bedge_type, const int *bedgeNum,
                            const double *x, const double *y, double *bc) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 7;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 7;
  }

  if(*bedge_type == 0) {
    for(int i = 0; i < 7; i++) {
      double y1 = y[exInd + i];
      // dBC[exInd + i] += y1 * (1.0 - y1);
      bc[exInd + i] += 2.0 * y1 * y1 * y1  - 3.0 * y1 * y1 + 1.0;
    }
  }
}
