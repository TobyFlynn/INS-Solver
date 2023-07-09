inline void mp_ins_2d_pr_1(const DG_FP *b0, const DG_FP *b1, const DG_FP *dt,
                           const DG_FP *dPdN, DG_FP *dPdNOld,
                           const DG_FP *divVelT, DG_FP *rho, DG_FP *pr_factor) {
  DG_FP factor = 1.0 / (*dt);
  for(int i = 0; i < DG_NP; i++) {
    divVelT[i] = -divVelT[i] * factor;
    pr_factor[i] = 1.0 / rho[i];
  }

  for(int i = 0; i < DG_NUM_FACES * DG_NPF; i++) {
    dPdNOld[i] = (*b0) * dPdN[i] + (*b1) * dPdNOld[i];
  }
}
