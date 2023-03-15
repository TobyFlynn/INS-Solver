inline void ins_advec_bc_2d(const DG_FP *t, const int *p, const int *bedge_type,
                            const int *bedgeNum, const DG_FP *x, const DG_FP *y,
                            const DG_FP *nx, const DG_FP *ny, const DG_FP *sJ,
                            const DG_FP *q0, const DG_FP *q1, DG_FP *flux0,
                            DG_FP *flux1) {
  const int exInd = *bedgeNum * DG_GF_NP;
  DG_FP pQ0[DG_GF_NP], pQ1[DG_GF_NP];
  const DG_FP PI = 3.141592653589793238463;

  // Get constants
  const DG_FP *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];

  // Set boundary velocities
  if(*bedge_type == 1) {
    // Outflow - Natural boundary condition
    for(int i = 0; i < DG_GF_NP; i++) {
      pQ0[i] = q0[exInd + i];
      pQ1[i] = q1[exInd + i];
    }
  } else {
    // Wall - slip
    if(fabs(y[exInd] - y[exInd + DG_GF_NP - 1]) < 1e-8) {
      for(int i = 0; i < DG_GF_NP; i++) {
        pQ0[i] = q0[exInd + i];
        pQ1[i] = 0.0;
      }
    } else {
      for(int i = 0; i < DG_GF_NP; i++) {
        pQ0[i] = 0.0;
        pQ1[i] = q1[exInd + i];
      }
    }
    // for(int i = 0; i < DG_GF_NP; i++) {
    //   const DG_FP dot = nx[exInd + i] * q0[exInd + i] + ny[exInd + i] * q1[exInd + i];
    //   pQ0[i] = q0[exInd + i] - dot * nx[exInd + i];
    //   pQ1[i] = q1[exInd + i] - dot * ny[exInd + i];
    // }
  }

  // Calculate numerical flux
  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_GF_NP; i++) {
    int ind = exInd + i;

    DG_FP mVel = q0[ind] * nx[ind] + q1[ind] * ny[ind];
    DG_FP pVel = pQ0[i] * -nx[ind] + pQ1[i] * -ny[ind];
    DG_FP vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Numerical flux calculation
  for(int i = 0; i < DG_GF_NP; i++) {
    int ind = exInd + i;
    // Get interior flux terms
    DG_FP mF0 = q0[ind] * q0[ind];
    DG_FP mF1 = q0[ind] * q1[ind];
    DG_FP mF2 = q0[ind] * q1[ind];
    DG_FP mF3 = q1[ind] * q1[ind];
    // Get exterior flux terms
    DG_FP pF0 = pQ0[i] * pQ0[i];
    DG_FP pF1 = pQ0[i] * pQ1[i];
    DG_FP pF2 = pQ0[i] * pQ1[i];
    DG_FP pF3 = pQ1[i] * pQ1[i];
    // Numerical flux
    flux0[ind] += 0.5 * gaussW[i] * sJ[ind] * (-nx[ind] * (mF0 - pF0) - ny[ind] * (mF1 - pF1) - maxVel * (pQ0[i] - q0[ind]));
    flux1[ind] += 0.5 * gaussW[i] * sJ[ind] * (-nx[ind] * (mF2 - pF2) - ny[ind] * (mF3 - pF3) - maxVel * (pQ1[i] - q1[ind]));
  }
}
