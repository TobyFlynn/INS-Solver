inline void cg3_update_res(const DG_FP *alpha, DG_FP *residual, const DG_FP *Ap0, 
                           const DG_FP *Ap1, const DG_FP *Ap2, DG_FP *r0, 
                           DG_FP *r1, DG_FP *r2) {
  DG_FP res_tmp = 0.0;
  for(int i = 0; i < DG_NP; i++) {
    r0[i] -= *alpha * Ap0[i];
    r1[i] -= *alpha * Ap1[i];
    r2[i] -= *alpha * Ap2[i];
    res_tmp += r0[i] * r0[i] + r1[i] * r1[i] + r2[i] * r2[i];
  }
  *residual += res_tmp;
}