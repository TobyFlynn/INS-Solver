inline void ls_advec_bedges(const int *p, const int *bedge_type,
                            const int *bedgeNum, const double *x,
                            const double *y, const double *q, double *exQ) {
  const int exInd = *bedgeNum * DG_GF_NP;

  for(int i = 0; i < DG_GF_NP; i++) {
    exQ[exInd + i] += q[exInd + i];
  }
}
