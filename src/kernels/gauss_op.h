inline void gauss_op(const double *tau, const double *sJ,
                     const double *mD0, double *f0_0, double *f0_1, double *f0_2,
                     const double *mD1, double *f1_0, double *f1_1, double *f1_2,
                     const double *mD2, double *f2_0, double *f2_1, double *f2_2,
                     double *pDy0, double *pDy1, double *pDy2) {
  // Face 0 - Transpose
  for(int ind = 0; ind < DG_GF_NP * DG_NP; ind++) {
    int indT = ((ind * DG_NP) % (DG_NP * DG_GF_NP)) + (ind / DG_GF_NP);
    f0_0[ind] = gFInterp0_g[indT];
    f0_1[ind] = gFInterp0_g[indT];
    f0_2[ind] = mD0[indT];
  }

  // Post diagonal multiply
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_GF_NP; n++) {
      int ind  = m * DG_GF_NP + n;
      f0_0[ind] = gaussW_g[n] * sJ[n] * tau[0] * f0_0[ind];
      f0_1[ind] = gaussW_g[n] * sJ[n] * f0_1[ind];
      f0_2[ind] = gaussW_g[n] * sJ[n] * f0_2[ind];
    }
  }

  // Face 1 - Transpose
  for(int ind = 0; ind < DG_GF_NP * DG_NP; ind++) {
    int indT = ((ind * DG_NP) % (DG_NP * 7)) + (ind / DG_GF_NP);
    f1_0[ind] = gFInterp1_g[indT];
    f1_1[ind] = gFInterp1_g[indT];
    f1_2[ind] = mD1[indT];
  }

  // Post diagonal multiply
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_GF_NP; n++) {
      int ind = m * DG_GF_NP + n;
      f1_0[ind] = gaussW_g[n] * sJ[n + DG_GF_NP] * tau[1] * f1_0[ind];
      f1_1[ind] = gaussW_g[n] * sJ[n + DG_GF_NP] * f1_1[ind];
      f1_2[ind] = gaussW_g[n] * sJ[n + DG_GF_NP] * f1_2[ind];
    }
  }

  // Face 2 - Transpose
  for(int ind = 0; ind < DG_GF_NP * DG_NP; ind++) {
    int indT = ((ind * DG_NP) % (DG_NP * DG_GF_NP)) + (ind / DG_GF_NP);
    f2_0[ind] = gFInterp2_g[indT];
    f2_1[ind] = gFInterp2_g[indT];
    f2_2[ind] = mD2[indT];
  }

  // Post diagonal multiply
  for(int m = 0; m < DG_NP; m++) {
    for(int n = 0; n < DG_GF_NP; n++) {
      int ind = m * DG_GF_NP + n;
      f2_0[ind] = gaussW_g[n] * sJ[n + 2 * DG_GF_NP] * tau[2] * f2_0[ind];
      f2_1[ind] = gaussW_g[n] * sJ[n + 2 * DG_GF_NP] * f2_1[ind];
      f2_2[ind] = gaussW_g[n] * sJ[n + 2 * DG_GF_NP] * f2_2[ind];
    }
  }

  for(int i = 0; i < DG_GF_NP * DG_NP; i++) {
    pDy0[i] = 0.0;
    pDy1[i] = 0.0;
    pDy2[i] = 0.0;
  }
}
