inline void art_vis_2d_max(const DG_FP *mu0, DG_FP *mu1) {
  for(int i = 0; i < DG_NP; i++) {
    if(mu0[i] > mu1[i])
      mu1[i] = mu0[i];
  }
}