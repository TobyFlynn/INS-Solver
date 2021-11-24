inline void ls_advec_bedges(const int *p, const int *bedge_type,
                            const int *bedgeNum, const double *x,
                            const double *y, const double *q, double *exQ) {
  // Get constants for this element's order
  const int dg_npf = DG_CONSTANTS[(*p - 1) * 5 + 1];
  const int *fmask = &FMASK[(*p - 1) * 3 * DG_NPF];
  fmask = &fmask[*bedgeNum * dg_npf];
  const int exInd = *bedgeNum * dg_npf;

  for(int i = 0; i < dg_npf; i++) {
    exQ[exInd + i] += q[fmask[i]];
  }
}
