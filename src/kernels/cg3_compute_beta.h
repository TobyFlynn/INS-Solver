inline void cg3_compute_beta(DG_FP *beta, const DG_FP *r0, const DG_FP *r1, const DG_FP *r2) {
  DG_FP tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    tmp += r0[i] * r0[i] + r1[i] * r1[i] + r2[i] * r2[i];
  }
  *beta += tmp;
}