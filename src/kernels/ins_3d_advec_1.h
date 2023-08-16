inline void ins_3d_advec_1(const int *faceNum, const int *fmaskL_corrected,
                           const int *fmaskR_corrected, const DG_FP *nx,
                           const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP **u,
                           const DG_FP **v, const DG_FP **w, DG_FP **f0,
                           DG_FP **f1, DG_FP **f2) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskL = &fmask[faceNum[0] * DG_NPF];
  const int *fmaskR = &fmask[faceNum[1] * DG_NPF];
  const int fIndL = faceNum[0] * DG_NPF;
  const int fIndR = faceNum[1] * DG_NPF;

  // Left numerical flux calculation
  const DG_FP int_factL = 0.5 * fscale[0];
  const DG_FP nxL = nx[0];
  const DG_FP nyL = ny[0];
  const DG_FP nzL = nz[0];
  const DG_FP nxR = nx[1];
  const DG_FP nyR = ny[1];
  const DG_FP nzR = nz[1];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskL_ind = fmaskL[i];
    const int fmaskR_ind = fmaskR_corrected[i];
    const DG_FP uM = u[0][fmaskL_ind];
    const DG_FP vM = v[0][fmaskL_ind];
    const DG_FP wM = w[0][fmaskL_ind];
    const DG_FP uP = u[1][fmaskR_ind];
    const DG_FP vP = v[1][fmaskR_ind];
    const DG_FP wP = w[1][fmaskR_ind];

    const DG_FP lVel = nxL * uM + nyL * vM + nzL * wM;
    const DG_FP rVel = nxR * uP + nyR * vP + nzR * wP;
    const DG_FP maxVel = fmax(fabs(lVel), fabs(rVel));

    f0[0][fIndL + i] = int_factL * (-nxL * (uM * uM - uP * uP)
                        - nyL * (uM * vM - uP * vP) - nzL * (uM * wM - uP * wP)
                        - maxVel * (uP - uM));
    f1[0][fIndL + i] = int_factL * (-nxL * (vM * uM - vP * uP)
                        - nyL * (vM * vM - vP * vP) - nzL * (vM * wM - vP * wP)
                        - maxVel * (vP - vM));
    f2[0][fIndL + i] = int_factL * (-nxL * (wM * uM - wP * uP)
                        - nyL * (wM * vM - wP * vP) - nzL * (wM * wM - wP * wP)
                        - maxVel * (wP - wM));
  }

  // Right numerical flux calculation
  const DG_FP int_factR = 0.5 * fscale[1];
  for(int i = 0; i < DG_NPF; i++) {
    const int fmaskR_ind = fmaskR[i];
    const int fmaskL_ind = fmaskL_corrected[i];
    const DG_FP uM = u[1][fmaskR_ind];
    const DG_FP vM = v[1][fmaskR_ind];
    const DG_FP wM = w[1][fmaskR_ind];
    const DG_FP uP = u[0][fmaskL_ind];
    const DG_FP vP = v[0][fmaskL_ind];
    const DG_FP wP = w[0][fmaskL_ind];

    const DG_FP lVel = nxL * uP + nyL * vP + nzL * wP;
    const DG_FP rVel = nxR * uM + nyR * vM + nzR * wM;
    const DG_FP maxVel = fmax(fabs(lVel), fabs(rVel));

    f0[1][fIndR + i] = int_factR * (-nxR * (uM * uM - uP * uP)
                        - nyR * (uM * vM - uP * vP) - nzR * (uM * wM - uP * wP)
                        - maxVel * (uP - uM));
    f1[1][fIndR + i] = int_factR * (-nxR * (vM * uM - vP * uP)
                        - nyR * (vM * vM - vP * vP) - nzR * (vM * wM - vP * wP)
                        - maxVel * (vP - vM));
    f2[1][fIndR + i] = int_factR * (-nxR * (wM * uM - wP * uP)
                        - nyR * (wM * vM - wP * vP) - nzR * (wM * wM - wP * wP)
                        - maxVel * (wP - wM));
  }
}
