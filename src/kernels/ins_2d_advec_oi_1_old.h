inline void ins_2d_advec_oi_1_old(const int *faceNum, const bool *reverse,
                  const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale, 
                  const DG_FP **u, const DG_FP **v, DG_FP **fluxU,
                  DG_FP **fluxV) {
  const bool rev = *reverse;

  // Left numerical flux calculation
  const DG_FP int_factL = 0.5 * fscale[0];
  const DG_FP nxL = nx[0];
  const DG_FP nyL = ny[0];
  const DG_FP nxR = nx[1];
  const DG_FP nyR = ny[1];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const int indL = faceNum[0] * DG_CUB_SURF_2D_NP + i;
    const int indR = rev ? (faceNum[1] + 1) * DG_CUB_SURF_2D_NP - i - 1 : faceNum[1] * DG_CUB_SURF_2D_NP + i;
    
    const DG_FP uM = u[0][indL];
    const DG_FP vM = v[0][indL];
    const DG_FP uP = u[1][indR];
    const DG_FP vP = v[1][indR];

    const DG_FP lVel = nxL * uM + nyL * vM;
    const DG_FP rVel = nxR * uP + nyR * vP;
    const DG_FP maxVel = fmax(fabs(lVel), fabs(rVel));

    fluxU[0][indL] = int_factL * (nxL * (uM * uM + uP * uP)
                   + nyL * (uM * vM + uP * vP) + maxVel * (uM - uP));
    fluxV[0][indL] = int_factL * (nxL * (vM * uM + vP * uP)
                   + nyL * (vM * vM + vP * vP) + maxVel * (vM - vP));
  }

  // Right numerical flux calculation
  const DG_FP int_factR = 0.5 * fscale[1];
  for(int i = 0; i < DG_CUB_SURF_2D_NP; i++) {
    const int indR = faceNum[1] * DG_CUB_SURF_2D_NP + i;
    const int indL = rev ? (faceNum[0] + 1) * DG_CUB_SURF_2D_NP - i - 1 : faceNum[0] * DG_CUB_SURF_2D_NP + i;
    
    const DG_FP uM = u[1][indR];
    const DG_FP vM = v[1][indR];
    const DG_FP uP = u[0][indL];
    const DG_FP vP = v[0][indL];

    const DG_FP lVel = nxL * uP + nyL * vP;
    const DG_FP rVel = nxR * uM + nyR * vM;
    const DG_FP maxVel = fmax(fabs(lVel), fabs(rVel));

    fluxU[1][indR] = int_factR * (nxR * (uM * uM + uP * uP)
                   + nyR * (uM * vM + uP * vP) + maxVel * (uM - uP));
    fluxV[1][indR] = int_factR * (nxR * (vM * uM + vP * uP)
                   + nyR * (vM * vM + vP * vP) + maxVel * (vM - vP));
  }
}
