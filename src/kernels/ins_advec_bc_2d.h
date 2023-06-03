inline void ins_advec_bc_2d(const DG_FP *t, const int *bedge_type, const int *bedgeNum,
                            const DG_FP *nx, const DG_FP *ny, const DG_FP *fscale,
                            const DG_FP *x, const DG_FP *y, const DG_FP *u,
                            const DG_FP *v, DG_FP *flux0, DG_FP *flux1) {
  DG_FP pU[DG_NPF], pV[DG_NPF];
  const DG_FP PI = 3.141592653589793238463;
  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF + bedgeNum[0] * DG_NPF];

  // Set boundary velocities
  if(*bedge_type == 0) {
    // Inflow - BC function dependant on time
    for(int i = 0; i < DG_NPF; i++) {
      pU[i] = 1.0;
      pV[i] = 0.0;
    }
  } else if(*bedge_type == 1) {
    // Outflow - Natural boundary condition
    for(int i = 0; i < DG_NPF; i++) {
      const int fmask_ind = fmask[i];
      pU[i] = u[fmask_ind];
      pV[i] = v[fmask_ind];
    }
  } else {
    // Wall - No slip
    for(int i = 0; i < DG_NPF; i++) {
      pU[i] = 0.0;
      pV[i] = 0.0;
    }
  }

  // Calculate numerical flux
  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];

    DG_FP mVel = u[fmask_ind] * nx[0] + v[fmask_ind] * ny[0];
    DG_FP pVel = pU[i] * -nx[0] + pV[i] * -ny[0];
    DG_FP vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Numerical flux calculation
  for(int i = 0; i < DG_NPF; i++) {
    const int fmask_ind = fmask[i];
    // Get interior flux terms
    DG_FP mF0 = u[fmask_ind] * u[fmask_ind];
    DG_FP mF1 = u[fmask_ind] * v[fmask_ind];
    DG_FP mF2 = u[fmask_ind] * v[fmask_ind];
    DG_FP mF3 = v[fmask_ind] * v[fmask_ind];
    // Get exterior flux terms
    DG_FP pF0 = pU[i] * pU[i];
    DG_FP pF1 = pU[i] * pV[i];
    DG_FP pF2 = pU[i] * pV[i];
    DG_FP pF3 = pV[i] * pV[i];
    // Numerical flux
    const int flux_ind = bedgeNum[0] * DG_NPF + i;
    flux0[flux_ind] = 0.5 * fscale[0] * (-nx[0] * (mF0 - pF0) - ny[0] * (mF1 - pF1) - maxVel * (pU[i] - u[fmask_ind]));
    flux1[flux_ind] = 0.5 * fscale[0] * (-nx[0] * (mF2 - pF2) - ny[0] * (mF3 - pF3) - maxVel * (pV[i] - v[fmask_ind]));
  }
}
