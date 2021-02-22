inline void pressure_bc3(double *pRHSex, double *pRHS) {
  for(int i = 0; i < 15; i++) {
    pRHS[FMASK[i]] += pRHSex[i];
    pRHSex[i] = 0.0;
  }
}
