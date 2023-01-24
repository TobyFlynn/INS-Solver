inline void euler_flux(const double *q0, const double *q1, 
                       const double *q2, const double *q3,
                       double *f0, double *f1, double *f2, double *f3,
                       double *g0, double *g1, double *g2, double *g3) {
  for(int i = 0; i < DG_NP; i++) {
    const double u = q1[i] / q0[i];
    const double v = q2[i] / q0[i];
    const double p = (gamma_e - 1.0) * (q3[i] - 0.5 * (q1[i] * u + q2[i] * v));

    f0[i] = q1[i];
    f1[i] = q1[i] * u + p;
    f2[i] = q2[i] * u;
    f3[i] = u * (q3[i] + p);

    g0[i] = q2[i];
    g1[i] = q1[i] * v;
    g2[i] = q2[i] * v + p;
    g3[i] = v * (q3[i] + p);
  }
}