inline void cg3_recp_factor_inv_mass(const DG_FP *factor, DG_FP *u, DG_FP *v, DG_FP *w) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] /= factor[i];
    v[i] /= factor[i];
    w[i] /= factor[i];
  }
}