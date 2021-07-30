inline void ls_advec_bedges(const int *bedge_type, const int *bedgeNum,
                            const double *x, const double *y, const double *q,
                            double *exQ) {
  int exInd = *bedgeNum * 5;
  int *fmask = &FMASK[*bedgeNum * 5];

  for(int i = 0; i < 5; i++) {
    exQ[exInd + i] += q[fmask[i]];
  }
}
