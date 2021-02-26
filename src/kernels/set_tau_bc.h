inline void set_tau_bc(const int *bedgeNum, const double *J, const double *sJ,
                       double *tau) {
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

  for(int i = 0; i < 5; i++) {
    tau[exInd + i] += 100.0*2.0*25.0*25.0 / (J[fmask[i]] / sJ[exInd + i]);
  }
}
