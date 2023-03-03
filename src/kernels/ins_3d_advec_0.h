inline void ins_3d_advec_0(const DG_FP *u, const DG_FP *v, const DG_FP *w,
                           DG_FP *f00, DG_FP *f01, DG_FP *f02, DG_FP *f10,
                           DG_FP *f11, DG_FP *f12, DG_FP *f20, DG_FP *f21,
                           DG_FP *f22) {
  for(int i = 0; i < DG_NP; i++) {
    f00[i] = u[i] * u[i]; f01[i] = u[i] * v[i]; f02[i] = u[i] * w[i];
    f10[i] = v[i] * u[i]; f11[i] = v[i] * v[i]; f12[i] = v[i] * w[i];
    f20[i] = w[i] * u[i]; f21[i] = w[i] * v[i]; f22[i] = w[i] * w[i];
  }
}
