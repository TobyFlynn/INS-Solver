inline void ls_2d_stencil_avg_1(const DG_FP *stencil, const DG_FP *avg, DG_FP *s) {
  if(stencil[0] != 1.0)
    return;

  const int *fmask = &FMASK[(DG_ORDER - 1) * DG_NUM_FACES * DG_NPF];
  for(int i = 0; i < DG_NUM_FACES * DG_NPF; i++) {
    if(avg[i] != 0.0) {
      int ind = fmask[i];
      s[i] = avg[i];
    }
  }
}