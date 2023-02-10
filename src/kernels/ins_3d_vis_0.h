inline void ins_3d_vis_0(const double *factor, double *u, double *v, double *w) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] *= *factor;
    v[i] *= *factor;
    w[i] *= *factor;
  }
}
