inline void poisson_test_bc(const int *bedge_type, const int *bedgeNum,
                            const double *x, const double *y, double *dBC) {
  int exInd = 0;
  if(*bedgeNum == 1) {
    exInd = 5;
  } else if(*bedgeNum == 2) {
    exInd = 2 * 5;
  }

  int *fmask;

  if(*bedgeNum == 0) {
    fmask = FMASK;
  } else if(*bedgeNum == 1) {
    fmask = &FMASK[5];
  } else {
    fmask = &FMASK[2 * 5];
  }

  if(*bedge_type == 0) {
    for(int i = 0; i < 5; i++) {
      double y1 = y[fmask[i]];
      // dBC[exInd + i] += y1 * (1.0 - y1);
      // dBC[exInd + i] += 2.0 * y1 * y1 * y1  - 3.0 * y1 * y1 + 1.0;
    }
  }
}
