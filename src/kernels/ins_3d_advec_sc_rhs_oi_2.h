inline void ins_3d_advec_sc_rhs_oi_2(const DG_FP *t, const int *bc_type, const int *faceNum,
                           const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *x, const DG_FP *y, const DG_FP *z,
                           const DG_FP *us, const DG_FP *vs, const DG_FP *ws,
                           const DG_FP *ub, const DG_FP *vb, const DG_FP *wb,
                           DG_FP *mUs, DG_FP *mVs, DG_FP *mWs,  DG_FP *pUs,
                           DG_FP *pVs, DG_FP *pWs, DG_FP *mUb, DG_FP *mVb,
                           DG_FP *mWb, DG_FP *pUb, DG_FP *pVb, DG_FP *pWb) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    mUs[fInd + i] = us[fmask_ind];
    mVs[fInd + i] = vs[fmask_ind];
    mWs[fInd + i] = ws[fmask_ind];
    mUb[fInd + i] = ub[fmask_ind];
    mVb[fInd + i] = vb[fmask_ind];
    mWb[fInd + i] = wb[fmask_ind];
  }

  if(*bc_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      pUs[fInd + i] = 0.0;
      pVs[fInd + i] = 0.0;
      pWs[fInd + i] = 0.0;
      pUb[fInd + i] = 0.0;
      pVb[fInd + i] = 0.0;
      pWb[fInd + i] = 0.0;
    }
  } else if(*bc_type == BC_TYPE_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      // S
      const DG_FP mags = sqrt(us[fmask_ind] * us[fmask_ind] + vs[fmask_ind] * vs[fmask_ind] + ws[fmask_ind] * ws[fmask_ind]);
      const DG_FP dots = *nx * us[fmask_ind] + *ny * vs[fmask_ind] + *nz * ws[fmask_ind];
      pUs[fInd + i] = us[fmask_ind] - dots * *nx;
      pVs[fInd + i] = vs[fmask_ind] - dots * *ny;
      pWs[fInd + i] = ws[fmask_ind] - dots * *nz;
      const DG_FP mag2s = sqrt(pUs[fInd + i] * pUs[fInd + i] + pVs[fInd + i] * pVs[fInd + i] + pWs[fInd + i] * pWs[fInd + i]);
      const DG_FP mag_factors = fabs(mags) < 1e-8 || fabs(mag2s) < 1e-8 ? 1.0 : mags / mag2s;
      pUs[fInd + i] *= mag_factors;
      pVs[fInd + i] *= mag_factors;
      pWs[fInd + i] *= mag_factors;
      // B
      const DG_FP magb = sqrt(ub[fmask_ind] * ub[fmask_ind] + vb[fmask_ind] * vb[fmask_ind] + wb[fmask_ind] * wb[fmask_ind]);
      const DG_FP dotb = *nx * ub[fmask_ind] + *ny * vb[fmask_ind] + *nz * wb[fmask_ind];
      pUb[fInd + i] = ub[fmask_ind] - dotb * *nx;
      pVb[fInd + i] = vb[fmask_ind] - dotb * *ny;
      pWb[fInd + i] = wb[fmask_ind] - dotb * *nz;
      const DG_FP mag2b = sqrt(pUb[fInd + i] * pUb[fInd + i] + pVb[fInd + i] * pVb[fInd + i] + pWb[fInd + i] * pWb[fInd + i]);
      const DG_FP mag_factorb = fabs(magb) < 1e-8 || fabs(mag2b) < 1e-8 ? 1.0 : magb / mag2b;
      pUb[fInd + i] *= mag_factorb;
      pVb[fInd + i] *= mag_factorb;
      pWb[fInd + i] *= mag_factorb;
    }
  } else if(*bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      pUs[fInd + i] = us[fmask_ind];
      pVs[fInd + i] = vs[fmask_ind];
      pWs[fInd + i] = ws[fmask_ind];
      pUb[fInd + i] = ub[fmask_ind];
      pVb[fInd + i] = vb[fmask_ind];
      pWb[fInd + i] = wb[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                             *nx, *ny, *nz, us[fmask_ind], vs[fmask_ind], ws[fmask_ind],
                             pUs[fInd + i], pVs[fInd + i], pWs[fInd + i]);
      ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                             *nx, *ny, *nz, ub[fmask_ind], vb[fmask_ind], wb[fmask_ind],
                             pUb[fInd + i], pVb[fInd + i], pWb[fInd + i]);
    }
  }
}
