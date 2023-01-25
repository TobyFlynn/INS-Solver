inline void euler_ic(const double *x, const double *y, double *q0,
                     double *q1, double *q2, double *q3) {
  const double PI = 3.141592653589793238463;
  const double R = 1.5;
  const double S = 13.5;
  const double M = 0.4;
  for(int i = 0; i < DG_NP; i++) {
    double f = (1.0 - x[i] * x[i] - y[i] * y[i]) / (2.0 * R * R);
    double u = (S * y[i] * exp(f)) / (2.0 * PI * R);
    double v = (1.0 - ((S * x[i] * exp(f)) / (2.0 * PI * R)));
    double rho_tmp = 1.0 - ((S * S * M * M * (gamma_e - 1.0) * exp(2.0 * f)) / (8.0 * PI * PI));
    double rho = pow(rho_tmp, 1.0 / (gamma_e - 1.0));
    double p = pow(rho, gamma_e) / (gamma_e * M * M);

    q0[i] = rho;
    q1[i] = rho * u;
    q2[i] = rho * v;
    q3[i] = 0.5 * (rho * u * u + rho * v * v) + p / (gamma_e - 1.0);
  }
}
