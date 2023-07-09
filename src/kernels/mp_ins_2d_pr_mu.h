inline void mp_ins_2d_pr_mu(const DG_FP *mu, DG_FP *in0) {
  for(int i = 0; i < DG_NP; i++) {
    in0[i] *= mu[i];
  }
}
