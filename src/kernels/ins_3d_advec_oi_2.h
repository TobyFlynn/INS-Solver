inline void ins_3d_advec_oi_2(const DG_FP *t, const int *bc_type, const int *faceNum,
                           const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                           const DG_FP *z, const DG_FP *u, const DG_FP *v,
                           const DG_FP *w, DG_FP *mU, DG_FP *mV, DG_FP *mW, 
                           DG_FP *pU, DG_FP *pV, DG_FP *pW) {
  const int *fmask  = &FMASK[(DG_ORDER - 1) * 4 * DG_NPF];
  const int *fmaskB = &fmask[*faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;
  const DG_FP PI = 3.141592653589793238463;
  if(*bc_type == LW_INFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      // uR[i] = sin(PI * (*t));
      DG_FP tmp = y[fmask_ind] * y[fmask_ind] + z[fmask_ind] * z[fmask_ind] - LW_INLET_RADIUS * LW_INLET_RADIUS;
      tmp = fabs(tmp);
      tmp = fmin(1.0, tmp / 0.1);
      mU[fInd + i] = u[fmask_ind];
      mV[fInd + i] = v[fmask_ind];
      mW[fInd + i] = w[fmask_ind];
      pU[fInd + i] = tmp * 1.0;
      pV[fInd + i] = 0.0;
      pW[fInd + i] = 0.0;
    }
  } else if(*bc_type == LW_OUTFLOW_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      mU[fInd + i] = u[fmask_ind];
      mV[fInd + i] = v[fmask_ind];
      mW[fInd + i] = w[fmask_ind];
      pU[fInd + i] = u[fmask_ind];
      pV[fInd + i] = v[fmask_ind];
      pW[fInd + i] = w[fmask_ind];
    }
  } else if(*bc_type == LW_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      const DG_FP mag = sqrt(u[fmask_ind] * u[fmask_ind] + v[fmask_ind] * v[fmask_ind] + w[fmask_ind] * w[fmask_ind]);
      const DG_FP dot = nx[0] * u[fmask_ind] + ny[0] * v[fmask_ind] + nz[0] * w[fmask_ind];
      pU[i] = u[fmask_ind] - dot * nx[0];
      pV[i] = v[fmask_ind] - dot * ny[0];
      pW[i] = w[fmask_ind] - dot * nz[0];
      const DG_FP mag2 = sqrt(pU[i] * pU[i] + pV[i] * pV[i] + pW[i] * pW[i]);
      const DG_FP mag_factor = fabs(mag) < 1e-8 || fabs(mag2) < 1e-8 ? 1.0 : mag / mag2;
      pU[i] *= mag_factor;
      pV[i] *= mag_factor;
      pW[i] *= mag_factor;
      mU[fInd + i] = u[fmask_ind];
      mV[fInd + i] = v[fmask_ind];
      mW[fInd + i] = w[fmask_ind];
    }
  } else if(*bc_type == LW_NO_SLIP_WALL_BC) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmaskB[i];
      pU[i] = 0.0;
      pV[i] = 0.0;
      pW[i] = 0.0;
      mU[fInd + i] = u[fmask_ind];
      mV[fInd + i] = v[fmask_ind];
      mW[fInd + i] = w[fmask_ind];
    }
  }
}
