inline void ins_3d_set_ic(const double *x, const double *y, const double *z,
                          double *u0, double *v0, double *w0, double *u1,
                          double *v1, double *w1) {
  for(int i = 0; i < DG_NP; i++) {
    // u0[i] = 0.0;
    // v0[i] = 0.0;
    // w0[i] = 0.0;
    u0[i] = sin(x[i]) * cos(y[i]) * cos(z[i]);
    v0[i] = -cos(x[i]) * sin(y[i]) * cos(z[i]);
    w0[i] = 0.0;
    u1[i] = u0[i];
    v1[i] = v0[i];
    w1[i] = w0[i];
  }
}
