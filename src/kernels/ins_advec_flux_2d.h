inline void ins_advec_flux_2d(const DG_FP *u, const DG_FP *v, DG_FP *f0,
                              DG_FP *f1, DG_FP *f2, DG_FP *f3) {
  for(int i = 0; i < DG_NP; i++) {
    f0[i] = u[i] * u[i];
    f1[i] = u[i] * v[i];
    f2[i] = u[i] * v[i];
    f3[i] = v[i] * v[i];
  }
}
