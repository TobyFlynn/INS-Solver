inline void calc_surf_ten(const double *nx, const double *ny, const double *curv,
                          const double *delta, double *sf_x, double *sf_y) {
  for(int i = 0; i < DG_NP; i++) {
    sf_x[i] = delta[i] * curv[i] * nx[i] / weber;
    sf_y[i] = delta[i] * curv[i] * ny[i] / weber;
  }
}