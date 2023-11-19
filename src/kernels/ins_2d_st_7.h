inline void ins_2d_st_7(const DG_FP *alpha_, const DG_FP *ls, DG_FP *stx, DG_FP *sty) {
  const DG_FP alpha = *alpha_;
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    DG_FP step = 0.5 * tanh(PI * ls[i] / alpha);
    DG_FP rho = 0.5 * rho0 * (1.0 + step) + 0.5 * rho1 * (1.0 - step);
    stx[i] /= rho;
    sty[i] /= rho;
  }
}
