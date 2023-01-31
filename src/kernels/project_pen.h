inline void project_pen(const double *factor, const double *u, const double *v,
                        const double *h, double *pen) {
  double max_vel = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    double vel = u[i] * u[i] + v[i] * v[i];
    if(vel > max_vel)
      max_vel = vel;
  }
  max_vel = sqrt(max_vel);

  *pen = *factor * (*h) * max_vel;
}
