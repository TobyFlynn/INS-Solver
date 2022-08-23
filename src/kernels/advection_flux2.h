inline void advection_flux2(const double *u0, const double *v0, 
                            const double *u1, const double *v1, 
                            double *f0, double *f1, double *f2, 
                            double *f3) {
  for(int i = 0; i < DG_NP; i++) {
    f0[i] = u0[i] * u1[i];
    f1[i] = u0[i] * v1[i];
    f2[i] = v0[i] * u1[i];
    f3[i] = v0[i] * v1[i];
  }
}
