inline void ins_3d_advec_sc_rhs_oi_2(const int *bc_type, const int *faceNum,
                           const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *x, const DG_FP *y, const DG_FP *z, 
                           const DG_FP *us, const DG_FP *vs, const DG_FP *ws, 
                           const DG_FP *ub, const DG_FP *vb, const DG_FP *wb, 
                           DG_FP *mUs, DG_FP *mVs, DG_FP *mWs,  DG_FP *pUs, 
                           DG_FP *pVs, DG_FP *pWs, DG_FP *mUb, DG_FP *mVb, 
                           DG_FP *mWb, DG_FP *pUb, DG_FP *pVb, DG_FP *pWb) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const DG_FP PI = 3.141592653589793238463;
  if(*bc_type == LW_INFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      // uR[i] = sin(PI * (*t));
      const int fmask_ind = fmaskB[i];
      DG_FP tmp = y[fmask_ind] * y[fmask_ind] + z[fmask_ind] * z[fmask_ind] - LW_INLET_RADIUS * LW_INLET_RADIUS;
      tmp = fabs(tmp);
      tmp = fmin(1.0, tmp / 0.1);
      mUs[i] = us[fmask_ind];
      mVs[i] = vs[fmask_ind];
      mWs[i] = ws[fmask_ind];
      mUb[i] = ub[fmask_ind];
      mVb[i] = vb[fmask_ind];
      mWb[i] = wb[fmask_ind];

      pUs[i] = tmp * 1.0;
      pVs[i] = 0.0;
      pWs[i] = 0.0;
      pUb[i] = tmp * 1.0;
      pVb[i] = 0.0;
      pWb[i] = 0.0;
    }
  } else if(*bc_type == LW_OUTFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      mUs[i] = us[fmask_ind];
      mVs[i] = vs[fmask_ind];
      mWs[i] = ws[fmask_ind];
      mUb[i] = ub[fmask_ind];
      mVb[i] = vb[fmask_ind];
      mWb[i] = wb[fmask_ind];

      pUs[i] = us[fmask_ind];
      pVs[i] = vs[fmask_ind];
      pWs[i] = ws[fmask_ind];
      pUb[i] = ub[fmask_ind];
      pVb[i] = vb[fmask_ind];
      pWb[i] = wb[fmask_ind];
    }
  } else if(*bc_type == LW_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      const DG_FP mags = sqrt(us[fmask_ind] * us[fmask_ind] + vs[fmask_ind] * vs[fmask_ind] + ws[fmask_ind] * ws[fmask_ind]);
      const DG_FP dots = nx[0] * us[fmask_ind] + ny[0] * vs[fmask_ind] + nz[0] * ws[fmask_ind];
      pUs[i] = us[fmask_ind] - dots * nx[0];
      pVs[i] = vs[fmask_ind] - dots * ny[0];
      pWs[i] = ws[fmask_ind] - dots * nz[0];
      const DG_FP mag2s = sqrt(pUs[i] * pUs[i] + pVs[i] * pVs[i] + pWs[i] * pWs[i]);
      const DG_FP mag_factors = fabs(mags) < 1e-8 || fabs(mag2s) < 1e-8 ? 1.0 : mags / mag2s;
      pUs[i] *= mag_factors;
      pVs[i] *= mag_factors;
      pWs[i] *= mag_factors;
      const DG_FP magb = sqrt(ub[fmask_ind] * ub[fmask_ind] + vb[fmask_ind] * vb[fmask_ind] + wb[fmask_ind] * wb[fmask_ind]);
      const DG_FP dotb = nx[0] * ub[fmask_ind] + ny[0] * vb[fmask_ind] + nz[0] * wb[fmask_ind];
      pUb[i] = ub[fmask_ind] - dotb * nx[0];
      pVb[i] = vb[fmask_ind] - dotb * ny[0];
      pWb[i] = wb[fmask_ind] - dotb * nz[0];
      const DG_FP mag2b = sqrt(pUb[i] * pUb[i] + pVb[i] * pVb[i] + pWb[i] * pWb[i]);
      const DG_FP mag_factorb = fabs(magb) < 1e-8 || fabs(mag2b) < 1e-8 ? 1.0 : magb / mag2b;
      pUb[i] *= mag_factorb;
      pVb[i] *= mag_factorb;
      pWb[i] *= mag_factorb;

      mUs[i] = us[fmask_ind];
      mVs[i] = vs[fmask_ind];
      mWs[i] = ws[fmask_ind];
      mUb[i] = ub[fmask_ind];
      mVb[i] = vb[fmask_ind];
      mWb[i] = wb[fmask_ind];
    }
  } else if(*bc_type == LW_NO_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      mUs[i] = us[fmask_ind];
      mVs[i] = vs[fmask_ind];
      mWs[i] = ws[fmask_ind];
      mUb[i] = ub[fmask_ind];
      mVb[i] = vb[fmask_ind];
      mWb[i] = wb[fmask_ind];

      pUs[i] = 0.0;
      pVs[i] = 0.0;
      pWs[i] = 0.0;
      pUb[i] = 0.0;
      pVb[i] = 0.0;
      pWb[i] = 0.0;
    }
  }
}
