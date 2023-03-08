inline void ls_step(const DG_FP *alpha, const DG_FP *s, DG_FP *rho,
                    DG_FP *mu) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    DG_FP step = tanh(PI * s[i] / *alpha);
    // DG_FP step = 0.0;
    // if(s[i] > *alpha) {
    //   step = 1.0;
    // } else if(s[i] > -*alpha) {
    //   step = 0.5 * (1.0 + (s[i] / *alpha) + (1.0 / PI) * sin(PI * s[i] / *alpha));
    //   // step = 0.5 * (1.0 + (1.0 / PI) * sin(PI * s[i] / *alpha));
    // }
    rho[i]  = 0.5 * rho0 * (1.0 + step) + 0.5 * rho1 * (1.0 - step);
    mu[i]   = 0.5 * mu0 * (1.0 + step) + 0.5 * mu1 * (1.0 - step);
    // rho[i]  = rho0 * step + (rho1 / rho0) * (1.0 - step);
    // mu[i]   = mu0 * step + (mu1 / mu0) * (1.0 - step);
    // rho[i] = 1.0;
    // mu[i] = 1.0;
  }
}
