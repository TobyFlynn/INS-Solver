inline void ls_advec_bedges(const int *bedge_type, const int *bedgeNum,
                            const double *x, const double *y, const double *q,
                            double *exQ) {
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

  for(int i = 0; i < 5; i++) {
    exQ[exInd + i] += q[fmask[i]];
  }
}
