inline void ls_step(const double *alpha, const double *s, double *step,
                    double *nu, double *rho) {
  const double PI = 3.141592653589793238463;
  for(int i = 0; i < 15; i++) {
    step[i] = tanh(PI * s[i] / *alpha);
    nu[i] = nu0 * step[i] + nu1 * (1.0 - step[i]);
    rho[i] = rho0 * step[i] + rho1 * (1.0 - step[i]);
  }
}
