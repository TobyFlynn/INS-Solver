inline void ins_3d_vis_0(const DG_FP *factor, DG_FP *u, DG_FP *v, DG_FP *w) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] *= *factor;
    v[i] *= *factor;
    w[i] *= *factor;
  }
}
