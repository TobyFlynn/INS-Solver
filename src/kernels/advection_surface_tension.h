inline void advection_surface_tension(const double *curv, const double *nx,
                                      const double *ny, const double *step_s,
                                      double *st_x, double *st_y) {
  for(int i = 0; i < DG_NP; i++) {
    // Calculate smoothed delta function
    // (only applying surface tension on the interface)
    double delta = 1.0 - fabs(step_s[i]);
    if(fabs(step_s[i]) >= 1.0) {
      delta = 0.0;
    }
    st_x[i] = delta * (curv[i] * nx[i] / weber);
    st_y[i] = delta * (curv[i] * ny[i] / weber);
  }
}
