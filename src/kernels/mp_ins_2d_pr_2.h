inline void mp_ins_2d_pr_2(const DG_FP *rho, DG_FP *dpdx, DG_FP *dpdy) {
  for(int i = 0; i < DG_NP; i++) {
    dpdx[i] /= rho[i];
    dpdy[i] /= rho[i];
  }
}
