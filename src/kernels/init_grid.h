inline void init_grid(const double **nc, double *nodeX, double *nodeY,
                      const double *xr, const double *yr, const double *xs,
                      const double *ys, double *rx, double *ry, double *sx,
                      double *sy, double *nx, double *ny, double *fscale) {
  nodeX[0] = nc[0][0];
  nodeX[1] = nc[1][0];
  nodeX[2] = nc[2][0];
  nodeY[0] = nc[0][1];
  nodeY[1] = nc[1][1];
  nodeY[2] = nc[2][1];

  // J = -xs.*yr + xr.*ys
  double J[15];
  for(int i = 0; i < 15; i++) {
    J[i] = -xs[i] * yr[i] + xr[i] * ys[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < 15; i++) {
    rx[i] = ys[i] / J[i];
    sx[i] = -yr[i] / J[i];
    ry[i] = -xs[i] / J[i];
    sy[i] = xr[i] / J[i];
  }

  // Calculate normals

  // Face 0
  for(int i = 0; i < 5; i++) {
    nx[i] = yr[FMASK[i]];
    ny[i] = -xr[FMASK[i]];
  }
  // Face 1
  for(int i = 0; i < 5; i++) {
    nx[5 + i] = ys[FMASK[5 + i]] - yr[FMASK[5 + i]];
    ny[5 + i] = xr[FMASK[5 + i]] - xs[FMASK[5 + i]];
  }
  // Face 2
  for(int i = 0; i < 5; i++) {
    nx[2 * 5 + i] = -ys[FMASK[2 * 5 + i]];
    ny[2 * 5 + i] = xs[FMASK[2 * 5 + i]];
  }

  // Normalise
  for(int i = 0; i < 3 * 5; i++) {
    double sJ = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / sJ;
    ny[i] = ny[i] / sJ;
    fscale[i] = sJ / J[FMASK[i]];
  }
}
