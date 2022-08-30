inline void advection_bc2(const int *p, const int *bedge_type,
                          const int *bedgeNum, const double *t, const int *problem,
                          const double *x, const double *y, const double *u0,
                          const double *v0, const double *u1, const double *v1,
                          const double *nx, const double *ny, const double *sJ, 
                          double *flux0, double *flux1) {
  const int exInd = *bedgeNum * DG_GF_NP;
  double pQ0[DG_GF_NP], pQ1[DG_GF_NP];
  const double PI = 3.141592653589793238463;

  // Get constants
  const double *gaussW = &gaussW_g[(*p - 1) * DG_GF_NP];

  // Set boundary velocities
  if(*problem == 0) {
    if(*bedge_type == 0) {
      // Inflow - BC function dependant on time
      for(int i = 0; i < DG_GF_NP; i++) {
        double y1 = y[exInd + i];
        // pQ0[i] = sin((PI * *t) / 8.0) * 4.0 * y1 * (1.0 - y1);
        pQ0[i] = 0.1;
        pQ1[i] = 0.0;
      }
    } else if(*bedge_type == 1) {
      // Outflow - Natural boundary condition
      for(int i = 0; i < DG_GF_NP; i++) {
        pQ0[i] = u1[exInd + i];
        pQ1[i] = v1[exInd + i];
      }
    } else {
      // Wall - No slip
      for(int i = 0; i < DG_GF_NP; i++) {
        pQ0[i] = 0.0;
        pQ1[i] = 0.0;
      }
    }
  } else {
    if(*bedge_type == 0) {
      // Inflow - BC function dependant on time
      for(int i = 0; i < DG_GF_NP; i++) {
        double y1 = y[exInd + i];
        double x1 = x[exInd + i];
        pQ0[i] += -sin(2.0 * PI * y1) * exp(-nu * 4.0 * PI * PI * *t);
        pQ1[i] += sin(2.0 * PI * x1) * exp(-nu * 4.0 * PI * PI * *t);
      }
    }

    if(*bedge_type == 1) {
      // Outflow - BC function dependant on time
      for(int i = 0; i < DG_GF_NP; i++) {
        pQ0[i] += u1[exInd + i];
        pQ1[i] += v1[exInd + i];
      }
    }
  }

  // Calculate numerical flux
  // Get max vel across the face
  int maxVel = 0.0;
  for(int i = 0; i < DG_GF_NP; i++) {
    int ind = exInd + i;

    double mVel = u1[ind] * nx[ind] + v1[ind] * ny[ind];
    double pVel = pQ0[i] * -nx[ind] + pQ1[i] * -ny[ind];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > maxVel) maxVel = vel;
  }

  // Numerical flux calculation
  for(int i = 0; i < DG_GF_NP; i++) {
    int ind = exInd + i;
    // Get interior flux terms
    double mF0 = u0[ind] * u1[ind];
    double mF1 = u0[ind] * v1[ind];
    double mF2 = v0[ind] * u1[ind];
    double mF3 = v0[ind] * v1[ind];
    // Get exterior flux terms
    double pF0 = u0[ind] * pQ0[i];
    double pF1 = u0[ind] * pQ1[i];
    double pF2 = v0[ind] * pQ1[i];
    double pF3 = v0[ind] * pQ1[i];
    // Numerical flux
    flux0[ind] += 0.5 * gaussW[i] * sJ[ind] * (-nx[ind] * (mF0 - pF0) - ny[ind] * (mF1 - pF1) - maxVel * (pQ0[i] - u1[ind]));
    flux1[ind] += 0.5 * gaussW[i] * sJ[ind] * (-nx[ind] * (mF2 - pF2) - ny[ind] * (mF3 - pF3) - maxVel * (pQ1[i] - v1[ind]));
  }
}
