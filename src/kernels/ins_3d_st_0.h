inline void ins_3d_st_0(const DG_FP *alpha_, DG_FP *out0) {
  const DG_FP alpha = *alpha_;
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_CUB_3D_NP; i++) {
    DG_FP step = 0.5 * tanh(PI * out0[i] / alpha);
    out0[i] = step;
  }
}
