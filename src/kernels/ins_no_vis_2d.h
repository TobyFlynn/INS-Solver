inline void ins_no_vis_2d(const DG_FP *g0, const DG_FP *velTT0,
                          const DG_FP *velTT1, DG_FP *velX, DG_FP *velY) {
  for(int i = 0; i < DG_NP; i++) {
    velX[i] = velTT0[i] / *g0;
    velY[i] = velTT1[i] / *g0;
  }
}
