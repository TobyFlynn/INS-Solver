inline void pressure_rhs(const double *b0, const double *b1, const double *dt,
                         const int *p, const double *J, const double *sJ, const double *dPdN,
                         double *dPdNOld, double *divVelT) {
  double factor = 1.0 / (*dt);
  for(int i = 0; i < DG_NP; i++) {
    // divVelT[i] = -divVelT[i] * factor * J[i];
    divVelT[i] = -divVelT[i] * factor;
  }

  const double *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];
  for(int i = 0; i < DG_G_NP; i++) {
    dPdNOld[i] = gaussW[i % DG_GF_NP] * sJ[i] * ((*b0) * dPdN[i] + (*b1) * dPdNOld[i]);
  }
}
