inline void art_vis_2d_max(const DG_FP *m0, DG_FP *m1) {
  for(int i = 0; i < DG_NP; i++) {
    if(m0[i] > m1[i])
      m1[i] = m0[i];
  }
}