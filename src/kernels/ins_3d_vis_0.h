inline void ins_3d_vis_0(const DG_FP *factor, const DG_FP *u, const DG_FP *v, 
                         const DG_FP *w, DG_FP *rhs0, DG_FP *rhs1, DG_FP *rhs2) {
  for(int i = 0; i < DG_NP; i++) {
    rhs0[i] = u[i] * *factor;
    rhs1[i] = v[i] * *factor;
    rhs2[i] = w[i] * *factor;
  }
}
