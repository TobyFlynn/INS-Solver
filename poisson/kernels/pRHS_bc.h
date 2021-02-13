inline void pRHS_bc(const int *bedge_type, const int *bedgeNum, const double *U,
                    double *exU) {
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

  if(*bedge_type == 1) {
    // Outflow
    for(int i = 0; i < 5; i++) {
      exU[exInd + i] += -U[fmask[i]];
    }
  }
}
