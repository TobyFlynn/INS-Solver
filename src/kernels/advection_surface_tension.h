inline void advection_surface_tension(const double *delta, const double *curv,
                                      const double *nx, const double *ny,
                                      const double *s, double *st_x,
                                      double *st_y) {
  for(int i = 0; i < DG_NP; i++) {
    st_x[i] = delta[i] * (curv[i] * nx[i] / weber);
    st_y[i] = delta[i] * (curv[i] * ny[i] / weber);
  }
}
