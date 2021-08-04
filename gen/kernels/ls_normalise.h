inline void ls_normalise(double *nx, double *ny) {
  for(int i = 0; i < 10; i++) {
    double size = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] = nx[i] / size;
    ny[i] = ny[i] / size;
  }
}
