inline void mp_ins_2d_add_mu(const DG_FP *mu, const DG_FP *rho, DG_FP *art_vis) {
  for(int i = 0; i < DG_NP; i++) {
    art_vis[i] *= rho[i] * r_ynolds;
    art_vis[i] += mu[i];
  }
}
