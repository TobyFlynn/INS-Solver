inline void calc_h_3d_max(const DG_FP *fscale, DG_FP *h) {
  if(fscale[0] < *h)
    *h = fscale[0];
  if(fscale[1] < *h)
    *h = fscale[1];
}
