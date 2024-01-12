inline void mp_ins_pr_bc_dpdn_factors(const DG_FP *b0, const DG_FP *b1,
                              const DG_FP *dPdN, DG_FP *dPdNOld) {
  for(int i = 0; i < DG_NUM_FACES * DG_NPF; i++) {
    dPdNOld[i] = (*b0) * dPdN[i] + (*b1) * dPdNOld[i];
  }
}
