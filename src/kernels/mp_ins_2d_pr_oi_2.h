inline void mp_ins_2d_pr_oi_2(const DG_FP *alpha_, DG_FP *out) {
  const DG_FP alpha = *alpha_;
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_2D_NP; i++) {
    DG_FP step = 0.5 * tanh(PI * out[i] / alpha);
    DG_FP rho = 0.5 * rho0 * (1.0 + step) + 0.5 * rho1 * (1.0 - step);
    out[i] = 1.0 / rho;
  }
}
