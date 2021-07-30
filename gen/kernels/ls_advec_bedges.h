inline void ls_advec_bedges(const int *bedge_type, const int *bedgeNum,
                            const double *x, const double *y, const double *q,
                            double *exQ) {
  int exInd = *bedgeNum * 4;
  int *fmask = &FMASK[*bedgeNum * 4];

  for(int i = 0; i < 4; i++) {
    exQ[exInd + i] += q[fmask[i]];
  }
}
