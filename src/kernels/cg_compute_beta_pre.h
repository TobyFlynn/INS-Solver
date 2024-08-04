inline void cg_compute_beta_pre(DG_FP *beta, const DG_FP *r0, const DG_FP *r1, 
                                const DG_FP *z0, const DG_FP *z1) {
  DG_FP tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    tmp += r0[i] * z0[i] + r1[i] * z1[i];
  }
  *beta += tmp;
}