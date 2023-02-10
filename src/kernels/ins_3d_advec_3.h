inline void ins_3d_advec_3(const double *a0, const double *a1, const double *b0,
                           const double *b1, const double *dt, const double *u,
                           const double *v, const double *w, const double *u_old,
                           const double *v_old, const double *w_old,
                           const double *n0, const double *n1, const double *n2,
                           const double *n0_old, const double *n1_old,
                           const double *n2_old, double *uT, double *vT,
                           double *wT) {
  for(int i = 0; i < DG_NP; i++) {
    uT[i] = (*a0 * u[i] + *a1 * u_old[i]) - *dt * (*b0 * n0[i] + *b1 * n0_old[i]);
    vT[i] = (*a0 * v[i] + *a1 * v_old[i]) - *dt * (*b0 * n1[i] + *b1 * n1_old[i]);
    wT[i] = (*a0 * w[i] + *a1 * w_old[i]) - *dt * (*b0 * n2[i] + *b1 * n2_old[i]);
  }
}
