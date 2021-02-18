inline void viscosity_bc(const int *bedge_type, const int *bedgeNum,
                         double *vRHS0, double *vRHS1) {
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

  if(*bedge_type == 0) {
    // Inflow
    for(int i = 0; i < 5; i++) {
      vRHS0[exInd + i] += bc_u;
      vRHS0[exInd + i] += bc_v;
      // vRHS0[exInd + i] += 1.0;
      // vRHS1[exInd + i] += 1.0;
    }
  }
}
