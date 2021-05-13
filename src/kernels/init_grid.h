inline void init_grid(double *rx, double *ry, double *sx, double *sy,
                      double *nx, double *ny, double *J, double *sJ, double *fscale) {
  // Calculate normals
  // Face 0
  for(int i = 0; i < 5; i++) {
    nx[i] = ry[FMASK[i]];
    ny[i] = -rx[FMASK[i]];
  }
  // Face 1
  for(int i = 0; i < 5; i++) {
    nx[5 + i] = sy[FMASK[5 + i]] - ry[FMASK[5 + i]];
    ny[5 + i] = rx[FMASK[5 + i]] - sx[FMASK[5 + i]];
  }
  // Face 2
  for(int i = 0; i < 5; i++) {
    nx[2 * 5 + i] = -sy[FMASK[2 * 5 + i]];
    ny[2 * 5 + i] = sx[FMASK[2 * 5 + i]];
  }

  // J = -xs.*yr + xr.*ys
  for(int i = 0; i < 15; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < 15; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  // Normalise
  for(int i = 0; i < 3 * 5; i++) {
    sJ[i] = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / sJ[i];
    ny[i] = ny[i] / sJ[i];
    fscale[i] = sJ[i] / J[FMASK[i]];
  }
}
