inline void poisson_rhs_qbc(const int *bedge_type, const int *bedgeNum,
                            const int *neumann0, const int *neumann1, const double *q, double *exq) {
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

  if(*bedge_type == *neumann0 || *bedge_type == *neumann1) {
    for(int i = 0; i < 5; i++) {
      exq[exInd + i] += -q[fmask[i]];
    }
  }
}