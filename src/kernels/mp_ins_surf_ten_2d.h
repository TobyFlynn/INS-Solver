inline void mp_ins_surf_ten_2d(const double *nx, const double *ny,
                               const double *curv, const double *delta_x, const double *delta_y,
                               double *sf_x, double *sf_y) {
  for(int i = 0; i < DG_NP; i++) {
    // sf_x[i] = 0.5 * fabs(delta_x[i]) * curv[i] * nx[i] / weber;
    // sf_y[i] = 0.5 * fabs(delta_y[i]) * curv[i] * ny[i] / weber;
    sf_x[i] = 0.5 * delta_x[i] * curv[i] / weber;
    sf_y[i] = 0.5 * delta_y[i] * curv[i] / weber;
  }
}
