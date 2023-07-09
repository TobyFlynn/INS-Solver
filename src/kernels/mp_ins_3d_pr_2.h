inline void mp_ins_3d_pr_2(const DG_FP *rho, DG_FP *dpdx, DG_FP *dpdy,
                           DG_FP *dpdz) {
  for(int i = 0; i < DG_NP; i++) {
    dpdx[i] /= rho[i];
    dpdy[i] /= rho[i];
    dpdz[i] /= rho[i];
  }
}
