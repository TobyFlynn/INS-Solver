inline void ls_advec_bedges(const int *bedge_type, const int *bedgeNum,
                            const double *x, const double *y, const double *q,
                            double *exQ) {
  int exInd = *bedgeNum * DG_NPF;
  int *fmask = &FMASK[*bedgeNum * DG_NPF];

  for(int i = 0; i < DG_NPF; i++) {
    exQ[exInd + i] += q[fmask[i]];
  }
}
