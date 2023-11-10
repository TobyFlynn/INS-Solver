inline void ins_3d_advec_oi_2(const DG_FP *t, const int *bc_type, const int *faceNum,
                           const DG_FP *nx, const DG_FP *ny, const DG_FP *nz,
                           const DG_FP *fscale, const DG_FP *x, const DG_FP *y,
                           const DG_FP *z, const DG_FP *u, const DG_FP *v,
                           const DG_FP *w, DG_FP *mU, DG_FP *mV, DG_FP *mW,
                           DG_FP *pU, DG_FP *pV, DG_FP *pW) {
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + *faceNum * DG_NPF];
  const int fInd = *faceNum * DG_NPF;

  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    mU[fInd + i] = u[fmask_ind];
    mV[fInd + i] = v[fmask_ind];
    mW[fInd + i] = w[fmask_ind];
  }

  if(*bc_type == BC_TYPE_NO_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      pU[fInd + i] = 0.0;
      pV[fInd + i] = 0.0;
      pW[fInd + i] = 0.0;
    }
  } else if(*bc_type == BC_TYPE_SLIP) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      const DG_FP mag = sqrt(u[fmask_ind] * u[fmask_ind] + v[fmask_ind] * v[fmask_ind] + w[fmask_ind] * w[fmask_ind]);
      const DG_FP dot = *nx * u[fmask_ind] + *ny * v[fmask_ind] + *nz * w[fmask_ind];
      pU[fInd + i] = u[fmask_ind] - dot * *nx;
      pV[fInd + i] = v[fmask_ind] - dot * *ny;
      pW[fInd + i] = w[fmask_ind] - dot * *nz;
      const DG_FP mag2 = sqrt(pU[fInd + i] * pU[fInd + i] + pV[fInd + i] * pV[fInd + i] + pW[fInd + i] * pW[fInd + i]);
      const DG_FP mag_factor = fabs(mag) < 1e-8 || fabs(mag2) < 1e-8 ? 1.0 : mag / mag2;
      pU[fInd + i] *= mag_factor;
      pV[fInd + i] *= mag_factor;
      pW[fInd + i] *= mag_factor;
    }
  } else if(*bc_type == BC_TYPE_NATURAL_OUTFLOW) {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      pU[fInd + i] = u[fmask_ind];
      pV[fInd + i] = v[fmask_ind];
      pW[fInd + i] = w[fmask_ind];
    }
  } else {
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      ps3d_custom_bc_get_vel(*bc_type, *t, x[fmask_ind], y[fmask_ind], z[fmask_ind],
                             *nx, *ny, *nz, u[fmask_ind], v[fmask_ind], w[fmask_ind],
                             pU[fInd + i], pV[fInd + i], pW[fInd + i]);
    }
  }
}
