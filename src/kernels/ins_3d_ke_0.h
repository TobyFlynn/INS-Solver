inline void ins_3d_ke_0(const DG_FP *u, const DG_FP *v, const DG_FP *w, DG_FP *ke) {
  for(int i = 0; i < DG_NP; i++) {
    ke[i] = u[i] * u[i] + v[i] * v[i] + w[i] * w[i];
  }
}