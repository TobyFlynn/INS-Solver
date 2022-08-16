inline void normalise_normals(double *nx, double *ny) {
  for(int i = 0; i < DG_NP; i++) {
    double mag = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
    nx[i] /= mag;
    ny[i] /= mag;
  }
}