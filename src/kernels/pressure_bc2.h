inline void pressure_bc2(const int *bedge_type, const int *bedgeNum,
                         const double *pRHS, double *pRHSex) {
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
    // Outflow - p = 0
    for(int i = 0; i < 5; i++) {
      int fInd = fmask[i];
      pRHSex[exInd + i] += -pRHS[fInd];
    }
  }
}
