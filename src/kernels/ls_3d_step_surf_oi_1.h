inline void ls_3d_step_surf_oi_1(const DG_FP *alpha, DG_FP *rho) {
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_3D_NP; i++) {
    DG_FP step = tanh(PI * rho[i] / *alpha);
    rho[i] = 0.5 * rho0 * (1.0 + step) + 0.5 * rho1 * (1.0 - step);
  }
}
