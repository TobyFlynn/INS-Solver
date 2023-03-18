inline void mp_ins_3d_surf_ten(const double *nx, const double *ny, const double *nz,
                               const double *curv, const double *delta_x,
                               const double *delta_y, const double *delta_z,
                               double *sf_x, double *sf_y, double *sf_z) {
  for(int i = 0; i < DG_NP; i++) {
    // sf_x[i] = 0.5 * fabs(delta_x[i]) * curv[i] * nx[i] / weber;
    // sf_y[i] = 0.5 * fabs(delta_y[i]) * curv[i] * ny[i] / weber;
    sf_x[i] = 0.5 * delta_x[i] * curv[i] / weber;
    sf_y[i] = 0.5 * delta_y[i] * curv[i] / weber;
    sf_z[i] = 0.5 * delta_z[i] * curv[i] / weber;
  }
}
