inline void ins_3d_advec_0(const double *u, const double *v, const double *w,
                           double *f00, double *f01, double *f02, double *f10,
                           double *f11, double *f12, double *f20, double *f21,
                           double *f22) {
  for(int i = 0; i < DG_NP; i++) {
    f00[i] = u[i] * u[i]; f01[i] = u[i] * v[i]; f02[i] = u[i] * w[i];
    f10[i] = v[i] * u[i]; f11[i] = v[i] * v[i]; f12[i] = v[i] * w[i];
    f20[i] = w[i] * u[i]; f21[i] = w[i] * v[i]; f22[i] = w[i] * w[i];
  }
}
