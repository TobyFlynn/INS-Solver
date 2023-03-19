inline void mp_ins_3d_pr_mu(const DG_FP *mu, DG_FP *in0, DG_FP *in1, DG_FP *in2) {
  for(int i = 0; i < DG_NP; i++) {
    in0[i] *= mu[i];
    in1[i] *= mu[i];
    in2[i] *= mu[i];
  }
}
