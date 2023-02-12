inline void mp_ins_3d_vis_0(const double *factor, const double *factor2,
                            const double *rho, double *u, double *v, double *w,
                            double *mm_factor) {
  for(int i = 0; i < DG_NP; i++) {
    u[i] *= *factor * rho[i];
    v[i] *= *factor * rho[i];
    w[i] *= *factor * rho[i];
    mm_factor[i] = *factor2 * rho[i];
  }
}