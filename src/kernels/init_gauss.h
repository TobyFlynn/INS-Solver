inline void init_gauss(double *rx, double *sx, double *ry, double *sy,
                       double *nx, double *ny, double *sJ) {

  // J = -xs.*yr + xr.*ys
  double J[21];
  for(int i = 0; i < 21; i++) {
    J[i] = -sx[i] * ry[i] + rx[i] * sy[i];
  }

  // rx = ys./J; sx =-yr./J; ry =-xs./J; sy = xr./J;
  for(int i = 0; i < 21; i++) {
    double rx_n = sy[i] / J[i];
    double sx_n = -ry[i] / J[i];
    double ry_n = -sx[i] / J[i];
    double sy_n = rx[i] / J[i];
    rx[i] = rx_n;
    sx[i] = sx_n;
    ry[i] = ry_n;
    sy[i] = sy_n;
  }

  // Calculate normals
  // Face 0
  for(int i = 0; i < 7; i++) {
    nx[i] = -sx[i];
    ny[i] = -sy[i];
  }
  // Face 1
  for(int i = 7; i < 14; i++) {
    nx[i] = rx[i] + sx[i];
    ny[i] = ry[i] + sy[i];
  }
  // Face 2
  for(int i = 14; i < 21; i++) {
    nx[i] = -rx[i];
    ny[i] = -ry[i];
  }

  // Normalise
  for(int i = 0; i < 21; i++) {
    sJ[i] = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / sJ[i];
    ny[i] = ny[i] / sJ[i];
    sJ[i] = sJ[i] * J[i];
  }
}
