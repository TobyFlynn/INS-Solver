inline void ls_3d_step_rho_vol_oi(const DG_FP *alpha, DG_FP *rho) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_CUB_3D_NP; i++) {
    DG_FP step = tanh(PI * rho[i] / *alpha);
    rho[i] = 0.5 * rho0 * (1.0 + step) + 0.5 * rho1 * (1.0 - step);
  }
}
