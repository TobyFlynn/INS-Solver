inline void mp_ins_2d_pr_grav(const DG_FP *avg_rho, const DG_FP *rho,
                              DG_FP *dpdx, DG_FP *dpdy) {
  const DG_FP grav_term = -1.0 / (froude * froude);
  for(int i = 0; i < DG_NP; i++) {
    dpdy[i] -= grav_term * (rho[i] - *avg_rho);
  }
}
