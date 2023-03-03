inline void euler_2d_faces(const int **p, const int *edgeNum, const bool *rev,
                           const DG_FP **nx, const DG_FP **ny, const DG_FP **sJ,
                           const DG_FP **gQ0, const DG_FP **gQ1, const DG_FP **gQ2,
                           const DG_FP **gQ3, DG_FP **gRHS0, DG_FP **gRHS1,
                           DG_FP **gRHS2, DG_FP **gRHS3) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;

  // Get constants
  const DG_FP *gaussW = &gaussW_g[(p[0][0] - 1) * DG_GF_NP];

  const int exIndL = edgeL * DG_GF_NP;
  const int exIndR = edgeR * DG_GF_NP;

  int maxLambda = 0.0;
  for(int i = 0; i < DG_GF_NP; i++) {
    int lInd = exIndL + i;
    int rInd;
    if(reverse) {
      rInd = exIndR + DG_GF_NP - 1 - i;
    } else {
      rInd = exIndR + i;
    }
    const DG_FP uL = gQ1[0][lInd] / gQ0[0][lInd];
    const DG_FP vL = gQ2[0][lInd] / gQ0[0][lInd];
    const DG_FP pL = (gamma_e - 1.0) * (gQ3[0][lInd] - 0.5 * (gQ1[0][lInd] * uL + gQ2[0][lInd] * vL));
    const DG_FP uR = gQ1[1][rInd] / gQ0[1][rInd];
    const DG_FP vR = gQ2[1][rInd] / gQ0[1][rInd];
    const DG_FP pR = (gamma_e - 1.0) * (gQ3[1][rInd] - 0.5 * (gQ1[1][rInd] * uR + gQ2[1][rInd] * vR));

    DG_FP lLambda = sqrt(uL * uL + vL * vL) + sqrt(fabs(gamma_e * pL / gQ0[0][lInd]));
    DG_FP rLambda = sqrt(uR * uR + vR * vR) + sqrt(fabs(gamma_e * pR / gQ0[1][rInd]));
    DG_FP lambda = fmax(lLambda, rLambda);
    if(lambda > maxLambda) maxLambda = lambda;
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
    const DG_FP uL = gQ1[0][lInd] / gQ0[0][lInd];
    const DG_FP vL = gQ2[0][lInd] / gQ0[0][lInd];
    const DG_FP pL = (gamma_e - 1.0) * (gQ3[0][lInd] - 0.5 * (gQ1[0][lInd] * uL + gQ2[0][lInd] * vL));
    const DG_FP uR = gQ1[1][rInd] / gQ0[1][rInd];
    const DG_FP vR = gQ2[1][rInd] / gQ0[1][rInd];
    const DG_FP pR = (gamma_e - 1.0) * (gQ3[1][rInd] - 0.5 * (gQ1[1][rInd] * uR + gQ2[1][rInd] * vR));

    const DG_FP f0L = gQ1[0][lInd];
    const DG_FP f1L = gQ1[0][lInd] * uL + pL;
    const DG_FP f2L = gQ2[0][lInd] * uL;
    const DG_FP f3L = uL * (gQ3[0][lInd] + pL);
    const DG_FP g0L = gQ2[0][lInd];
    const DG_FP g1L = gQ1[0][lInd] * vL;
    const DG_FP g2L = gQ2[0][lInd] * vL + pL;
    const DG_FP g3L = vL * (gQ3[0][lInd] + pL);
    const DG_FP f0R = gQ1[1][rInd];
    const DG_FP f1R = gQ1[1][rInd] * uR + pR;
    const DG_FP f2R = gQ2[1][rInd] * uR;
    const DG_FP f3R = uR * (gQ3[1][rInd] + pR);
    const DG_FP g0R = gQ2[1][rInd];
    const DG_FP g1R = gQ1[1][rInd] * vR;
    const DG_FP g2R = gQ2[1][rInd] * vR + pR;
    const DG_FP g3R = vR * (gQ3[1][rInd] + pR);

    gRHS0[0][lInd] += 0.5 * gaussW[i] * sJ[0][lInd] * (nx[0][lInd] * (f0L + f0R) + ny[0][lInd] * (g0L + g0R) + maxLambda * (gQ0[0][lInd] - gQ0[1][rInd]));
    gRHS1[0][lInd] += 0.5 * gaussW[i] * sJ[0][lInd] * (nx[0][lInd] * (f1L + f1R) + ny[0][lInd] * (g1L + g1R) + maxLambda * (gQ1[0][lInd] - gQ1[1][rInd]));
    gRHS2[0][lInd] += 0.5 * gaussW[i] * sJ[0][lInd] * (nx[0][lInd] * (f2L + f2R) + ny[0][lInd] * (g2L + g2R) + maxLambda * (gQ2[0][lInd] - gQ2[1][rInd]));
    gRHS3[0][lInd] += 0.5 * gaussW[i] * sJ[0][lInd] * (nx[0][lInd] * (f3L + f3R) + ny[0][lInd] * (g3L + g3R) + maxLambda * (gQ3[0][lInd] - gQ3[1][rInd]));

    gRHS0[1][rInd] += 0.5 * gaussW[i] * sJ[1][rInd] * (nx[1][rInd] * (f0L + f0R) + ny[1][rInd] * (g0L + g0R) + maxLambda * (gQ0[1][rInd] - gQ0[0][lInd]));
    gRHS1[1][rInd] += 0.5 * gaussW[i] * sJ[1][rInd] * (nx[1][rInd] * (f1L + f1R) + ny[1][rInd] * (g1L + g1R) + maxLambda * (gQ1[1][rInd] - gQ1[0][lInd]));
    gRHS2[1][rInd] += 0.5 * gaussW[i] * sJ[1][rInd] * (nx[1][rInd] * (f2L + f2R) + ny[1][rInd] * (g2L + g2R) + maxLambda * (gQ2[1][rInd] - gQ2[0][lInd]));
    gRHS3[1][rInd] += 0.5 * gaussW[i] * sJ[1][rInd] * (nx[1][rInd] * (f3L + f3R) + ny[1][rInd] * (g3L + g3R) + maxLambda * (gQ3[1][rInd] - gQ3[0][lInd]));
  }
}