inline void ins_advec_faces_2d(const int **p, const int *edgeNum, const bool *rev,
                               const double **nx, const double **ny,
                               const double **sJ, const double **q0,
                               const double **q1, double **flux0,
                               double **flux1) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;

  // Get constants
  const double *gaussW = &gaussW_g[(p[0][0] - 1) * DG_GF_NP];

  const int exIndL = edgeL * DG_GF_NP;
  const int exIndR = edgeR * DG_GF_NP;

  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd = exIndL + i;
    int rInd;
    if(reverse) {
      rInd = exIndR + DG_GF_NP - 1 - i;
    } else {
      rInd = exIndR + i;
    }

    double lVel = q0[0][lInd] * nx[0][lInd] + q1[0][lInd] * ny[0][lInd];
    double rVel = q0[1][rInd] * nx[1][rInd] + q1[1][rInd] * ny[1][rInd];
    double vel = fmax(fabs(lVel), fabs(rVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Numerical flux calculation
  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd = exIndL + i;
    int rInd;
    if(reverse) {
      rInd = exIndR + DG_GF_NP - 1 - i;
    } else {
      rInd = exIndR + i;
    }
    // Get left flux terms
    double lF0 = q0[0][lInd] * q0[0][lInd];
    double lF1 = q0[0][lInd] * q1[0][lInd];
    double lF2 = q0[0][lInd] * q1[0][lInd];
    double lF3 = q1[0][lInd] * q1[0][lInd];
    // Get right flux terms
    double rF0 = q0[1][rInd] * q0[1][rInd];
    double rF1 = q0[1][rInd] * q1[1][rInd];
    double rF2 = q0[1][rInd] * q1[1][rInd];
    double rF3 = q1[1][rInd] * q1[1][rInd];
    // Numerical flux LHS
    flux0[0][lInd] += 0.5 * gaussW[i] * sJ[0][lInd] * (-nx[0][lInd] * (lF0 - rF0) - ny[0][lInd] * (lF1 - rF1) - maxVel * (q0[1][rInd] - q0[0][lInd]));
    flux1[0][lInd] += 0.5 * gaussW[i] * sJ[0][lInd] * (-nx[0][lInd] * (lF2 - rF2) - ny[0][lInd] * (lF3 - rF3) - maxVel * (q1[1][rInd] - q1[0][lInd]));
    // Numerical flux RHS
    flux0[1][rInd] += 0.5 * gaussW[i] * sJ[1][rInd] * (-nx[1][rInd] * (rF0 - lF0) - ny[1][rInd] * (rF1 - lF1) - maxVel * (q0[0][lInd] - q0[1][rInd]));
    flux1[1][rInd] += 0.5 * gaussW[i] * sJ[1][rInd] * (-nx[1][rInd] * (rF2 - lF2) - ny[1][rInd] * (rF3 - lF3) - maxVel * (q1[0][lInd] - q1[1][rInd]));
  }
}
