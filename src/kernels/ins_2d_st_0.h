inline void ins_2d_st_0(const DG_FP *alpha_, DG_FP *out0) {
  const DG_FP alpha = *alpha_;
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_CUB_2D_NP; i++) {
    // DG_FP step = 0.5 * tanh(PI * out0[i] / alpha);
    DG_FP step = 0.5 * (1 + out0[i] / alpha + sin(PI * out0[i] / alpha) / PI);
    if(out0[i] < -alpha) step = 0.0;
    if(out0[i] > alpha) step = 1.0;
    out0[i] = step;
  }
}