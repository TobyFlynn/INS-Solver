inline void ls_step_vol_oi(const DG_FP *_alpha, DG_FP *rho) {
  const DG_FP PI = 3.141592653589793238463;
  const DG_FP alpha = *_alpha;
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    // DG_FP step = tanh(PI * rho[i] / alpha);
    // rho[i] = 0.5 * rho0 * (1.0 + step) + 0.5 * rho1 * (1.0 - step);
    DG_FP step = 0.5 * (1 + rho[i] / alpha + sin(PI * rho[i] / alpha) / PI);
    if(rho[i] < -alpha) step = 0.0;
    if(rho[i] > alpha) step = 1.0;
    rho[i] = rho1 + (rho0 - rho1) * step;
  }
}