inline void set_ic2(const int *bedge_type, const int *bedgeNum,
                    double *uD, double *qN) {
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

  // Dirichlet BC u = 1
  if(*bedge_type == 0) {
    for(int i = 0; i < 5; i++) {
      // uD[exInd + i] += -1.0;
      uD[exInd + i] += 1.0;
    }
  } else if(*bedge_type == 1) {
    for(int i = 0; i < 5; i++) {
      // uD[exInd + i] += -0.75;
      uD[exInd + i] += 1.0;
    }
  } else if(*bedge_type == 2) {
    for(int i = 0; i < 5; i++) {
      // uD[exInd + i] += 0.5;
      uD[exInd + i] += 1.0;
    }
  } else {
    for(int i = 0; i < 5; i++) {
      // uD[exInd + i] += 1.2;
      uD[exInd + i] += 1.0;
    }
  }
}
