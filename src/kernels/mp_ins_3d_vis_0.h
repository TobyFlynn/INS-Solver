inline void mp_ins_3d_vis_0(const double *fact0, const double *fact1,
                            const double *rho, const double *mu, 
                            const double *art_vis, double *u, double *v, 
                            double *w, double *factor, double *mm_factor) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] *= *fact0 * rho[i];
    v[i] *= *fact0 * rho[i];
    w[i] *= *fact0 * rho[i];
    // factor[i] = mu[i] + *art_vis;
    factor[i] = mu[i];
    mm_factor[i] = *fact1 * rho[i];
  }
}