inline void poisson_mult_faces(const int *pL, const DG_FP *uL, const DG_FP *opL, DG_FP *rhsL,
                               const int *pR, const DG_FP *uR, const DG_FP *opR, DG_FP *rhsR) {
  const int dg_npL = DG_CONSTANTS[(*pL - 1) * DG_NUM_CONSTANTS];
  const int dg_npR = DG_CONSTANTS[(*pR - 1) * DG_NUM_CONSTANTS];

  for(int m = 0; m < dg_npL; m++) {
    for(int n = 0; n < dg_npR; n++) {
      int ind = m + n * dg_npL;
      rhsL[m] += opL[ind] * uR[n];
    }
  }

  for(int m = 0; m < dg_npR; m++) {
    for(int n = 0; n < dg_npL; n++) {
      int ind = m + n * dg_npR;
      rhsR[m] += opR[ind] * uL[n];
    }
  }
}
