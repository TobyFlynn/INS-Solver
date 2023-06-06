inline void ins_proj_setup_1(const DG_FP *fscale, DG_FP *h) {
  DG_FP max_fscale = 0.0;
  for(int i = 0; i < DG_NUM_FACES; i++) {
    max_fscale = fmax(max_fscale, fscale[i]);
  }
  *h = 1.0 / max_fscale;
}
