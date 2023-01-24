inline void euler_edges(const int **p, const int *edgeNum, const bool *rev,
                        const double **nx, const double **ny, const double **sJ,
                        const double **gQ0, const double **gQ1, const double **gQ2,
                        const double **gQ3, double **gRHS0, double **gRHS1, 
                        double **gRHS2, double **gRHS3) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;

  // Get constants
  const double *gaussW = &gaussW_g[(p[0][0] - 1) * DG_GF_NP];

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
    const double uL = gQ1[0][lInd] / gQ0[0][lInd];
    const double vL = gQ2[0][lInd] / gQ0[0][lInd];
    const double pL = (gamma_e - 1.0) * (gQ3[0][lInd] - 0.5 * (gQ1[0][lInd] * uL + gQ2[0][lInd] * vL));
    const double uR = gQ1[1][rInd] / gQ0[1][rInd];
    const double vR = gQ2[1][rInd] / gQ0[1][rInd];
    const double pR = (gamma_e - 1.0) * (gQ3[1][rInd] - 0.5 * (gQ1[1][rInd] * uR + gQ2[1][rInd] * vR));

    double lLambda = sqrt(uL * uL + vL * vL) + sqrt(fabs(gamma_e * pL / gQ0[0][lInd]));
    double rLambda = sqrt(uR * uR + vR * vR) + sqrt(fabs(gamma_e * pR / gQ0[1][rInd]));
    double lambda = fmax(lLambda, rLambda);
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
    const double uL = gQ1[0][lInd] / gQ0[0][lInd];
    const double vL = gQ2[0][lInd] / gQ0[0][lInd];
    const double pL = (gamma_e - 1.0) * (gQ3[0][lInd] - 0.5 * (gQ1[0][lInd] * uL + gQ2[0][lInd] * vL));
    const double uR = gQ1[1][rInd] / gQ0[1][rInd];
    const double vR = gQ2[1][rInd] / gQ0[1][rInd];
    const double pR = (gamma_e - 1.0) * (gQ3[1][rInd] - 0.5 * (gQ1[1][rInd] * uR + gQ2[1][rInd] * vR));

    const double f0L = gQ1[0][lInd];
    const double f1L = gQ1[0][lInd] * uL + pL;
    const double f2L = gQ2[0][lInd] * uL;
    const double f3L = uL * (gQ3[0][lInd] + pL);
    const double g0L = gQ2[0][lInd];
    const double g1L = gQ1[0][lInd] * vL;
    const double g2L = gQ2[0][lInd] * vL + pL;
    const double g3L = vL * (gQ3[0][lInd] + pL);
    const double f0R = gQ1[1][rInd];
    const double f1R = gQ1[1][rInd] * uR + pR;
    const double f2R = gQ2[1][rInd] * uR;
    const double f3R = uR * (gQ3[1][rInd] + pR);
    const double g0R = gQ2[1][rInd];
    const double g1R = gQ1[1][rInd] * vR;
    const double g2R = gQ2[1][rInd] * vR + pR;
    const double g3R = vR * (gQ3[1][rInd] + pR);

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