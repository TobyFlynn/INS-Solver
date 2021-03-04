inline void reset_nBC(double *nBCx, double *nBCy) {
  for(int i = 0; i < 15; i++) {
    nBCx[i] = 0.0;
    nBCy[i] = 0.0;
  }
}
