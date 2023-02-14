inline void mp_ins_3d_vis_0(const DG_FP *fact0, const DG_FP *fact1,
                            const DG_FP *rho, const DG_FP *mu, 
                            const DG_FP *art_vis, DG_FP *u, DG_FP *v, 
                            DG_FP *w, DG_FP *factor, DG_FP *mm_factor) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] *= *fact0 * rho[i];
    v[i] *= *fact0 * rho[i];
    w[i] *= *fact0 * rho[i];
    // factor[i] = mu[i] + *art_vis;
    factor[i] = mu[i];
    mm_factor[i] = *fact1 * rho[i];
  }
}