inline void ls_step(const double *alpha, const double *s, double *rho,
                    double *mu) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    double step = tanh(PI * s[i] / *alpha);
    // rho[i]  = 0.5 * rho0 * (1.0 + step[i]) + 0.5 * rho1 * (1.0 - step[i]);
    // mu[i]   = 0.5 * mu0 * (1.0 + step[i]) + 0.5 * mu1 * (1.0 - step[i]);
    rho[i] = 1.0;
    mu[i] = 1.0;
  }
}
