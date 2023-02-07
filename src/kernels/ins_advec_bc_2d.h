inline void ins_advec_bc_2d(const double *t, const int *p, const int *bedge_type,
                            const int *bedgeNum, const double *x, const double *y,
                            const double *nx, const double *ny, const double *sJ,
                            const double *q0, const double *q1, double *flux0,
                            double *flux1) {
  const int exInd = *bedgeNum * DG_GF_NP;
  double pQ0[DG_GF_NP], pQ1[DG_GF_NP];
  const double PI = 3.141592653589793238463;

  // Get constants
  const double *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];

  // Set boundary velocities
  if(*bedge_type == 0) {
    // Inflow - BC function dependant on time
    for(int i = 0; i < DG_GF_NP; i++) {
      pQ0[i] = 1.0;
      pQ1[i] = 0.0;
    }
  } else if(*bedge_type == 1) {
    // Outflow - Natural boundary condition
    for(int i = 0; i < DG_GF_NP; i++) {
      pQ0[i] = q0[exInd + i];
      pQ1[i] = q1[exInd + i];
    }
  } else {
    // Wall - No slip
    for(int i = 0; i < DG_GF_NP; i++) {
      pQ0[i] = 0.0;
      pQ1[i] = 0.0;
    }
  }

  // Calculate numerical flux
  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_GF_NP; i++) {
    int ind = exInd + i;

    double mVel = q0[ind] * nx[ind] + q1[ind] * ny[ind];
    double pVel = pQ0[i] * -nx[ind] + pQ1[i] * -ny[ind];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Numerical flux calculation
  for(int i = 0; i < DG_GF_NP; i++) {
    int ind = exInd + i;
    // Get interior flux terms
    double mF0 = q0[ind] * q0[ind];
    double mF1 = q0[ind] * q1[ind];
    double mF2 = q0[ind] * q1[ind];
    double mF3 = q1[ind] * q1[ind];
    // Get exterior flux terms
    double pF0 = pQ0[i] * pQ0[i];
    double pF1 = pQ0[i] * pQ1[i];
    double pF2 = pQ0[i] * pQ1[i];
    double pF3 = pQ1[i] * pQ1[i];
    // Numerical flux
    flux0[ind] += 0.5 * gaussW[i] * sJ[ind] * (-nx[ind] * (mF0 - pF0) - ny[ind] * (mF1 - pF1) - maxVel * (pQ0[i] - q0[ind]));
    flux1[ind] += 0.5 * gaussW[i] * sJ[ind] * (-nx[ind] * (mF2 - pF2) - ny[ind] * (mF3 - pF3) - maxVel * (pQ1[i] - q1[ind]));
  }
}