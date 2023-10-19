inline void ins_2d_st_art_vis(const DG_FP *max_diff_, const DG_FP *diff_width_,
                              const DG_FP *trans_width_, const DG_FP *alpha_,
                              const DG_FP *s, DG_FP *out) {
  const DG_FP alpha = *alpha_ * *diff_width_;
  const DG_FP max_diff = *max_diff_;
  const DG_FP trans_width = *alpha_ * *trans_width_;
  const DG_FP PI = 3.141592653589793238463;
  for(int i = 0; i < DG_NP; i++) {
    if(fabs(s[i]) <= alpha) {
      out[i] = max_diff;
    } else {
      out[i] = fmax(1e-6, fmin(max_diff, max_diff * (trans_width - (fabs(s[i]) - alpha)) / trans_width));
    }
  }
}
