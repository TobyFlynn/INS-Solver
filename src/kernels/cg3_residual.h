inline void cg3_residual(DG_FP *residual, const DG_FP *rhs0, const DG_FP *rhs1, 
                         const DG_FP *rhs2, DG_FP *r0,  DG_FP *r1, DG_FP *r2) {
  DG_FP res_tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    r0[i] = rhs0[i] - r0[i];
    r1[i] = rhs1[i] - r1[i];
    r2[i] = rhs2[i] - r2[i];
    res_tmp += r0[i] * r0[i] + r1[i] * r1[i] + r2[i] * r2[i];
  }
  *residual += res_tmp;
}