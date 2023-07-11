inline void euler_2d_faces(const int *edgeNum, const bool *rev,
                           const DG_FP *nx, const DG_FP *ny,
                           const DG_FP *fscale, const DG_FP **Q0,
                           const DG_FP **Q1, const DG_FP **Q2, const DG_FP **Q3,
                           DG_FP **RHS0, DG_FP **RHS1, DG_FP **RHS2,
                           DG_FP **RHS3) {
  // Work out which edge for each element
  const int edgeL = edgeNum[0];
  const int edgeR = edgeNum[1];
  const bool reverse = *rev;
  const int *fmaskL = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeL * DG_NPF];
  const int *fmaskR = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + edgeR * DG_NPF];

  int maxLambda = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = reverse ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    const DG_FP uL = Q1[0][fmaskL_ind] / Q0[0][fmaskL_ind];
    const DG_FP vL = Q2[0][fmaskL_ind] / Q0[0][fmaskL_ind];
    const DG_FP pL = (gamma_e - 1.0) * (Q3[0][fmaskL_ind] - 0.5 * (Q1[0][fmaskL_ind] * uL + Q2[0][fmaskL_ind] * vL));
    const DG_FP uR = Q1[1][fmaskR_ind] / Q0[1][fmaskR_ind];
    const DG_FP vR = Q2[1][fmaskR_ind] / Q0[1][fmaskR_ind];
    const DG_FP pR = (gamma_e - 1.0) * (Q3[1][fmaskR_ind] - 0.5 * (Q1[1][fmaskR_ind] * uR + Q2[1][fmaskR_ind] * vR));

    DG_FP lLambda = sqrt(uL * uL + vL * vL) + sqrt(fabs(gamma_e * pL / Q0[0][fmaskL_ind]));
    DG_FP rLambda = sqrt(uR * uR + vR * vR) + sqrt(fabs(gamma_e * pR / Q0[1][fmaskR_ind]));
    DG_FP lambda = fmax(lLambda, rLambda);
    if(lambda > maxLambda) maxLambda = lambda;
  }

  // Numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = reverse ? fmaskR[DG_NPF - i - 1] : fmaskR[i];

    const DG_FP uL = Q1[0][fmaskL_ind] / Q0[0][fmaskL_ind];
    const DG_FP vL = Q2[0][fmaskL_ind] / Q0[0][fmaskL_ind];
    const DG_FP pL = (gamma_e - 1.0) * (Q3[0][fmaskL_ind] - 0.5 * (Q1[0][fmaskL_ind] * uL + Q2[0][fmaskL_ind] * vL));
    const DG_FP uR = Q1[1][fmaskR_ind] / Q0[1][fmaskR_ind];
    const DG_FP vR = Q2[1][fmaskR_ind] / Q0[1][fmaskR_ind];
    const DG_FP pR = (gamma_e - 1.0) * (Q3[1][fmaskR_ind] - 0.5 * (Q1[1][fmaskR_ind] * uR + Q2[1][fmaskR_ind] * vR));

    const DG_FP f0L = Q1[0][fmaskL_ind];
    const DG_FP f1L = Q1[0][fmaskL_ind] * uL + pL;
    const DG_FP f2L = Q2[0][fmaskL_ind] * uL;
    const DG_FP f3L = uL * (Q3[0][fmaskL_ind] + pL);
    const DG_FP g0L = Q2[0][fmaskL_ind];
    const DG_FP g1L = Q1[0][fmaskL_ind] * vL;
    const DG_FP g2L = Q2[0][fmaskL_ind] * vL + pL;
    const DG_FP g3L = vL * (Q3[0][fmaskL_ind] + pL);
    const DG_FP f0R = Q1[1][fmaskR_ind];
    const DG_FP f1R = Q1[1][fmaskR_ind] * uR + pR;
    const DG_FP f2R = Q2[1][fmaskR_ind] * uR;
    const DG_FP f3R = uR * (Q3[1][fmaskR_ind] + pR);
    const DG_FP g0R = Q2[1][fmaskR_ind];
    const DG_FP g1R = Q1[1][fmaskR_ind] * vR;
    const DG_FP g2R = Q2[1][fmaskR_ind] * vR + pR;
    const DG_FP g3R = vR * (Q3[1][fmaskR_ind] + pR);

    const int outL_ind = edgeL * DG_NPF + i;
    const int outR_ind = reverse ? edgeR * DG_NPF + DG_NPF - i - 1: edgeR * DG_NPF + i;

    RHS0[0][outL_ind] = 0.5 * fscale[0] * (nx[0] * (f0L + f0R) + ny[0] * (g0L + g0R) + maxLambda * (Q0[0][fmaskL_ind] - Q0[1][fmaskR_ind]));
    RHS1[0][outL_ind] = 0.5 * fscale[0] * (nx[0] * (f1L + f1R) + ny[0] * (g1L + g1R) + maxLambda * (Q1[0][fmaskL_ind] - Q1[1][fmaskR_ind]));
    RHS2[0][outL_ind] = 0.5 * fscale[0] * (nx[0] * (f2L + f2R) + ny[0] * (g2L + g2R) + maxLambda * (Q2[0][fmaskL_ind] - Q2[1][fmaskR_ind]));
    RHS3[0][outL_ind] = 0.5 * fscale[0] * (nx[0] * (f3L + f3R) + ny[0] * (g3L + g3R) + maxLambda * (Q3[0][fmaskL_ind] - Q3[1][fmaskR_ind]));

    RHS0[1][outR_ind] = 0.5 * fscale[1] * (nx[1] * (f0L + f0R) + ny[1] * (g0L + g0R) + maxLambda * (Q0[1][fmaskR_ind] - Q0[0][fmaskL_ind]));
    RHS1[1][outR_ind] = 0.5 * fscale[1] * (nx[1] * (f1L + f1R) + ny[1] * (g1L + g1R) + maxLambda * (Q1[1][fmaskR_ind] - Q1[0][fmaskL_ind]));
    RHS2[1][outR_ind] = 0.5 * fscale[1] * (nx[1] * (f2L + f2R) + ny[1] * (g2L + g2R) + maxLambda * (Q2[1][fmaskR_ind] - Q2[0][fmaskL_ind]));
    RHS3[1][outR_ind] = 0.5 * fscale[1] * (nx[1] * (f3L + f3R) + ny[1] * (g3L + g3R) + maxLambda * (Q3[1][fmaskR_ind] - Q3[0][fmaskL_ind]));
  }
}
