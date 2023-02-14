inline void mp_ins_3d_pr_1(const DG_FP *b0, const DG_FP *b1, const DG_FP *dt,
                           const DG_FP *dPdN, DG_FP *dPdN_old, DG_FP *divVelT,
                           const DG_FP *rho, DG_FP *pr_factor) {
  DG_FP factor = 1.0 / *dt;
  for(int i = 0; i < DG_NP; i++) {
    divVelT[i] = -divVelT[i] * factor;
    pr_factor[i] = 1.0 / rho[i];
  }

  for(int i = 0; i < 4 * DG_NPF; i++) {
    dPdN_old[i] = (*b0) * dPdN[i] + (*b1) * dPdN_old[i];
  }
}