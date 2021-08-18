inline void advection_surface_tension(const double *alpha, const double *curv,
                                      const double *nx, const double *ny,
                                      const double *s, double *st_x,
                                      double *st_y) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 10; i++) {
    // Calculate smoothed delta function
    // (only applying surface tension on the interface)
    // double delta = 1.0 - fabs(step_s[i]);
    // double delta = 0.25 * (1.0 + cos(PI * step_s[i] / 4.0));
    // if(fabs(step_s[i]) >= 1.0) {
    //   delta = 0.0;
    // }
    // double delta = fabs(fabs(step_s[i]) - 1.0) / 2.0 * log(2.0);
    double delta = (PI / *alpha) * (1.0 / cosh(PI * s[i] / *alpha)) * (1.0 / cosh(PI * s[i] / *alpha));
    st_x[i] = delta * (curv[i] * nx[i] / weber);
    st_y[i] = delta * (curv[i] * ny[i] / weber);
  }
}
