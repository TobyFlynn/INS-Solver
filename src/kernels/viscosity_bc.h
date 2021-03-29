// TODO double check these
inline void viscosity_bc(const int *bedge_type, const int *bedgeNum,
                         const double *t, const double *x, const double *y,
                         double *exQ0, double *exQ1) {
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
    // Inflow - BC function dependant on time
    const double PI = 3.141592653589793238463;
    for(int i = 0; i < 5; i++) {
      double y1 = y[fmask[i]];
      exQ0[exInd + i] += pow(0.41, -2.0) * sin((PI * (*t)) / 8.0) * 6.0 * y1 * (0.41 - y1);
      // exQ0[exInd + i] += pow(0.41, -2.0) * cos((PI * *t) / 8.0);
      // exQ0[exInd + i] += 0.1;
      // exQ1[exInd + i] += bc_v;
    }
  }
}
