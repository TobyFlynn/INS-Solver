inline void mp_ins_2d_pr_bc_oi_1(const DG_FP *b0, const DG_FP *b1, const DG_FP *dt,
                           const DG_FP *dPdN, DG_FP *dPdNOld, DG_FP *divVelT) {
  DG_FP factor = 1.0 / (*dt);
  for(int i = 0; i < DG_NP; i++) {
    divVelT[i] = -divVelT[i] * factor;
  }

  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_2D_NP; i++) {
    dPdNOld[i] = (*b0) * dPdN[i] + (*b1) * dPdNOld[i];
  }
}
