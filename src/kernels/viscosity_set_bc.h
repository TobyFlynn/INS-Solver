inline void viscosity_set_bc(double *qtt0, double *qtt1, double *exQ0, double *exQ1) {
  for(int i = 0; i < 15; i++) {
    qtt0[FMASK[i]] = exQ0[i];
    qtt1[FMASK[i]] = exQ1[i];
    exQ0[i] = 0.0;
    exQ1[i] = 0.0;
  }
}
