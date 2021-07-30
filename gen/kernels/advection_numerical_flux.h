inline void advection_numerical_flux(const double *fscale, const double *nx,
                                     const double *ny, const double *q0,
                                     const double *q1, double *exQ0,
                                     double *exQ1, double *flux0, double *flux1) {
  // Compute fluxes for face nodes
  double fM[4][3 * 4];
  for(int i = 0; i < 3 * 4; i++) {
    fM[0][i] = q0[FMASK[i]] * q0[FMASK[i]];
    fM[1][i] = q0[FMASK[i]] * q1[FMASK[i]];
    fM[2][i] = q0[FMASK[i]] * q1[FMASK[i]];
    fM[3][i] = q1[FMASK[i]] * q1[FMASK[i]];
  }
  double fP[4][3 * 4];
  for(int i = 0; i < 3 * 4; i++) {
    fP[0][i] = exQ0[i] * exQ0[i];
    fP[1][i] = exQ0[i] * exQ1[i];
    fP[2][i] = exQ0[i] * exQ1[i];
    fP[3][i] = exQ1[i] * exQ1[i];
  }

  // Compute max velocity across each face
  double maxVel[3 * 4];
  double max = 0.0;
  for(int i = 0; i < 4; i++) {
    double mVel = q0[FMASK[i]] * nx[i] + q1[FMASK[i]] * ny[i];
    double pVel = exQ0[i] * nx[i] + exQ1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 0; i < 4; i++) {
    maxVel[i] = max;
  }
  max = 0.0;
  for(int i = 4; i < 2 * 4; i++) {
    double mVel = q0[FMASK[i]] * nx[i] + q1[FMASK[i]] * ny[i];
    double pVel = exQ0[i] * nx[i] + exQ1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 4; i < 2 * 4; i++) {
    maxVel[i] = max;
  }
  max = 0.0;
  for(int i = 2 * 4; i < 3 * 4; i++) {
    double mVel = q0[FMASK[i]] * nx[i] + q1[FMASK[i]] * ny[i];
    double pVel = exQ0[i] * nx[i] + exQ1[i] * ny[i];
    double vel = fmax(fabs(mVel), fabs(pVel));
    if(vel > max) max = vel;
  }
  for(int i = 2 * 4; i < 3 * 4; i++) {
    maxVel[i] = max;
  }

  // Lax-Friedrichs
  for(int i = 0; i < 3 * 4; i++) {
    flux0[i] = 0.5 * fscale[i] * (-nx[i] * (fM[0][i] - fP[0][i]) - ny[i] * (fM[1][i] - fP[1][i]) - maxVel[i] * (exQ0[i] - q0[FMASK[i]]));
    flux1[i] = 0.5 * fscale[i] * (-nx[i] * (fM[2][i] - fP[2][i]) - ny[i] * (fM[3][i] - fP[3][i]) - maxVel[i] * (exQ1[i] - q1[FMASK[i]]));
  }

  // Zero exQ
  for(int i = 0; i < 3 * 4; i++) {
    exQ0[i] = 0.0;
    exQ1[i] = 0.0;
  }
}
