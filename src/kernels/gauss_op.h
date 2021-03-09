inline void gauss_op(const double *tau, const double *sJ,
                     const double *mD0, double *f0_0, double *f0_1, double *f0_2,
                     const double *mD1, double *f1_0, double *f1_1, double *f1_2,
                     const double *mD2, double *f2_0, double *f2_1, double *f2_2) {
  // // Face 0
  // for(int m = 0; m < 7; m++) {
  //   for(int n = 0; n < 15; n++) {
  //     int ind = m * 15 + n;
  //     f0_0[ind] = gaussW[m] * sJ[m] * tau[0] * gFInterp0[ind];
  //     f0_1[ind] = gaussW[m] * sJ[m] * gFInterp0[ind];
  //     f0_2[ind] = gaussW[m] * sJ[m] * mD0[ind];
  //   }
  // }
  //
  // // Face 1
  // for(int m = 0; m < 7; m++) {
  //   for(int n = 0; n < 15; n++) {
  //     int ind = m * 15 + n;
  //     f1_0[ind] = gaussW[m] * sJ[m + 7] * tau[1] * gFInterp1[ind];
  //     f1_1[ind] = gaussW[m] * sJ[m + 7] * gFInterp1[ind];
  //     f1_2[ind] = gaussW[m] * sJ[m + 7] * mD1[ind];
  //   }
  // }
  //
  // // Face 2
  // for(int m = 0; m < 7; m++) {
  //   for(int n = 0; n < 15; n++) {
  //     int ind = m * 15 + n;
  //     f2_0[ind] = gaussW[m] * sJ[m + 14] * tau[2] * gFInterp2[ind];
  //     f2_1[ind] = gaussW[m] * sJ[m + 14] * gFInterp2[ind];
  //     f2_2[ind] = gaussW[m] * sJ[m + 14] * mD2[ind];
  //   }
  // }

  // Face 0 - Transpose
  for(int ind = 0; ind < 7 * 15; ind++) {
    int indT = (ind * 15) % (15 * 7) + (ind / 7);
    f0_0[ind] = gFInterp0[indT];
    f0_1[ind] = gFInterp0[indT];
    f0_2[ind] = mD0[indT];
  }

  // Post diagonal multiply
  for(int m = 0; m < 15; m++) {
    for(int n = 0; n < 7; n++) {
      int ind  = m * 7 + n;
      f0_0[ind] = gaussW[n] * sJ[n] * tau[0] * f0_0[ind];
      f0_1[ind] = gaussW[n] * sJ[n] * f0_1[ind];
      f0_2[ind] = gaussW[n] * sJ[n] * f0_2[ind];
    }
  }

  // Face 1 - Transpose
  for(int ind = 0; ind < 7 * 15; ind++) {
    int indT = (ind * 15) % (15 * 7) + (ind / 7);
    f1_0[ind] = gFInterp1[indT];
    f1_1[ind] = gFInterp1[indT];
    f1_2[ind] = mD1[indT];
  }

  // Post diagonal multiply
  for(int m = 0; m < 15; m++) {
    for(int n = 0; n < 7; n++) {
      int ind = m * 7 + n;
      f1_0[ind] = gaussW[n] * sJ[n + 7] * tau[1] * f1_0[ind];
      f1_1[ind] = gaussW[n] * sJ[n + 7] * f1_1[ind];
      f1_2[ind] = gaussW[n] * sJ[n + 7] * f1_2[ind];
    }
  }

  // Face 2 - Transpose
  for(int ind = 0; ind < 7 * 15; ind++) {
    int indT = (ind * 15) % (15 * 7) + (ind / 7);
    f2_0[ind] = gFInterp2[indT];
    f2_1[ind] = gFInterp2[indT];
    f2_2[ind] = mD2[indT];
  }

  // Post diagonal multiply
  for(int m = 0; m < 15; m++) {
    for(int n = 0; n < 7; n++) {
      int ind = m * 7 + n;
      f2_0[ind] = gaussW[n] * sJ[n + 14] * tau[2] * f2_0[ind];
      f2_1[ind] = gaussW[n] * sJ[n + 14] * f2_1[ind];
      f2_2[ind] = gaussW[n] * sJ[n + 14] * f2_2[ind];
    }
  }
}
