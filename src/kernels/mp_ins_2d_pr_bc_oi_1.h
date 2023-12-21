inline void mp_ins_2d_pr_bc_oi_1(const DG_FP *b0, const DG_FP *b1,
                                 const DG_FP *dPdN, DG_FP *dPdNOld) {
  for(int i = 0; i < DG_NUM_FACES * DG_CUB_SURF_2D_NP; i++) {
    dPdNOld[i] = (*b0) * dPdN[i] + (*b1) * dPdNOld[i];
  }
}
